"""Verified account sessions and request-scoped ownership; no client identity headers."""
from __future__ import annotations
from contextlib import contextmanager
from contextvars import ContextVar
import json
import os
import urllib.error
import urllib.request

from fastapi import HTTPException
from starlette.concurrency import run_in_threadpool
from starlette.responses import JSONResponse

principal: ContextVar[str | None] = ContextVar('scientific_account', default=None)

@contextmanager
def account_scope(user_id):
    token = principal.set(user_id)
    try:
        yield
    finally:
        principal.reset(token)


def verify_session(cookie: str) -> str:
    if not cookie or len(cookie) > 8192:
        raise HTTPException(401, '请先登录')
    url = os.getenv('MARINE_AUTH_URL', '')
    if not url:
        raise HTTPException(503, '账号验证服务尚未配置')
    request = urllib.request.Request(url, headers={'Cookie': cookie, 'Accept': 'application/json',
                                                 'Host': os.getenv('MARINE_AUTH_HOST', 'lebb.chat')})
    try:
        client = urllib.request.build_opener(urllib.request.ProxyHandler({}))
        with client.open(request, timeout=5) as response:
            payload = json.loads(response.read(65536))
        if not isinstance(payload, dict):
            raise HTTPException(503, '账号验证服务返回了无效响应')
        data = payload.get('data') if payload.get('code') == 0 else None
        user_id = data.get('id') if isinstance(data, dict) else None
        if not isinstance(user_id, (str, int)) or isinstance(user_id, bool) or not str(user_id).strip():
            raise HTTPException(401, '登录状态无效')
        return str(user_id)
    except urllib.error.HTTPError as exc:
        raise HTTPException(401 if exc.code in (401,403) else 503, '登录状态无效' if exc.code in (401,403) else '账号验证服务暂不可用') from exc
    except (OSError, ValueError) as exc:
        raise HTTPException(503, '账号验证服务暂不可用') from exc


class AccountMiddleware:
    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope['type'] != 'http' or scope['path'] == '/health' or scope['path'].startswith('/internal/'):
            return await self.app(scope, receive, send)
        headers = dict(scope.get('headers', []))
        try:
            user_id = await run_in_threadpool(verify_session, headers.get(b'cookie', b'').decode('latin1'))
        except HTTPException as exc:
            return await JSONResponse({'detail': exc.detail}, status_code=exc.status_code)(scope, receive, send)
        with account_scope(user_id):
            await self.app(scope, receive, send)

# The API installs durable ownership hooks; CLI/scientific workers have no
# account context. Tool execution captures the submitting account explicitly.
analysis_check = None
analysis_grant = None

def require_analysis(kind, analysis_id):
    if principal.get() is not None and (analysis_check is None or not analysis_check(kind, analysis_id, principal.get())):
        raise ValueError('分析记录不存在或不属于当前账号')

def grant_analysis(kind, analysis_id):
    if principal.get() is not None:
        if analysis_grant is None:
            raise RuntimeError('analysis ownership store is unavailable')
        analysis_grant(kind, analysis_id, principal.get())

def scoped_tools(tools):
    import functools
    owner = principal.get()
    if owner is None:
        return tools
    from .progress import reporter
    progress_callback = reporter.get()
    def scoped(function):
        @functools.wraps(function)
        def invoke(*args, **kwargs):
            with account_scope(owner):
                token = reporter.set(progress_callback)
                try:
                    return function(*args, **kwargs)
                finally:
                    reporter.reset(token)
        return invoke
    return [t.model_copy(update={'func': scoped(t.func)}) if getattr(t, 'func', None) else t for t in tools]
