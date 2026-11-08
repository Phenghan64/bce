"""Bounded transport for the official ColabFold endpoint, with optional SSH egress.

Only the network route changes. TLS verification, queries and database search
settings remain the same; no single-sequence fallback is performed.
"""
from __future__ import annotations

import os
import time
import json
from pathlib import Path
from urllib.parse import urlsplit

import requests


class MsaTransportFailure(SystemExit):
    """Exit the upstream client's otherwise unbounded timeout retry loop."""

    def __init__(self, code: str, reason: str, retryable: bool = True):
        super().__init__(1)
        self.detail = {"error_code": code, "reason": reason, "retryable": retryable}


class MsaTransport:
    def __init__(self, original, *, proxy: str | None = None, sleep=time.sleep, clock=time.monotonic):
        self.original = original
        self.proxy = proxy if proxy is not None else os.environ.get("OPENFOLD_MSA_PROXY", "")
        if proxy is None and not self.proxy:
            config = Path.home() / '.config/marine-agent/msa-network.json'
            if config.is_file():
                self.proxy = json.loads(config.read_text()).get('fallback_proxy', '')
        self.sleep, self.clock = sleep, clock
        self.started = clock()
        self.route = "direct"
        self.attempts: list[dict] = []
        if self.proxy:
            parsed = urlsplit(self.proxy)
            if parsed.scheme != "socks5h" or parsed.hostname not in {"127.0.0.1", "localhost"} or not parsed.port:
                raise ValueError("OPENFOLD_MSA_PROXY must reference a loopback socks5h SSH tunnel")

    def request(self, session, method, url, **kwargs):
        parsed = urlsplit(url)
        if parsed.hostname != "api.colabfold.com" or parsed.scheme != "https":
            return self.original(session, method, url, **kwargs)
        if self.clock() - self.started > 600:
            raise MsaTransportFailure("msa_deadline", "MSA 服务等待超时，请稍后继续任务。")
        method = method.upper()
        for attempt in range(3):
            route = self.route
            options = {**kwargs, "timeout": (10, 120)}
            options["proxies"] = {"http": self.proxy, "https": self.proxy} if route == "ssh_egress" else {"http": None, "https": None, "all": None}
            row = {"method": method, "path": parsed.path, "route": route}
            self.attempts.append(row)
            print(f"MSA network: {method} {parsed.path} · {route} · attempt {attempt + 1}/3", flush=True)
            try:
                response = self.original(session, method, url, **options)
            except requests.exceptions.SSLError:
                row["error"] = "tls_verification"
                raise MsaTransportFailure("msa_tls", "MSA 服务 TLS 校验失败，请检查网络证书。", False)
            except requests.exceptions.ReadTimeout:
                row["error"] = "read_timeout"
                if method == "POST":
                    raise MsaTransportFailure("msa_submission_uncertain", "MSA 提交响应超时，未确认服务端是否受理；本次不自动重复提交。")
                response = None
            except (requests.exceptions.ConnectionError, requests.exceptions.ConnectTimeout):
                row["error"] = "connection"
                response = None
                if self.proxy:
                    self.route = "ssh_egress"
            if response is not None:
                row["http"] = response.status_code
                if response.status_code < 400:
                    return response
                if response.status_code not in {429, 502, 503, 504}:
                    raise MsaTransportFailure("msa_http", f"MSA 服务返回 HTTP {response.status_code}。", False)
                if method == "POST" and response.status_code != 429:
                    raise MsaTransportFailure("msa_submission_uncertain", f"MSA 提交返回 HTTP {response.status_code}，未自动重复提交。")
                try:
                    delay = min(15, max(1, float(response.headers.get("Retry-After", "1"))))
                except ValueError:
                    delay = 2
                response.close()
            else:
                delay = min(4, 2 ** attempt)
            if attempt < 2:
                self.sleep(delay)
        raise MsaTransportFailure("msa_network", "MSA 服务连接失败，已完成有限重试；可稍后继续，序列和输入已保留。")

    def report(self):
        return {"service": "https://api.colabfold.com", "route": self.route, "attempts": self.attempts,
                "tls_verified": True, "single_sequence_fallback": False}
