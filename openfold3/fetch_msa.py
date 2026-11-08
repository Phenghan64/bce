"""Fetch/cache an example query's MSA with a longer network read timeout."""
import json
import sys
from pathlib import Path

import requests
from msa_transport import MsaTransport, MsaTransportFailure

transport = MsaTransport(requests.sessions.Session.request)


def request_with_longer_timeout(self, method, url, **kwargs):
    return transport.request(self, method, url, **kwargs)


requests.sessions.Session.request = request_with_longer_timeout

from openfold3.core.data.tools.colabfold_msa_server import (  # noqa: E402
    MsaComputationSettings,
    preprocess_colabfold_msas,
)
from openfold3.projects.of3_all_atom.config.inference_query_format import (  # noqa: E402
    InferenceQuerySet,
)

query_path = Path(sys.argv[1]).resolve()
output = Path(sys.argv[2]).resolve()
output.mkdir(parents=True, exist_ok=True)
settings = MsaComputationSettings.from_config_with_cli_override(output, None)
settings._set_workspace_root(output / "openfold3_intermediates" / "colabfold_msas")
queries = InferenceQuerySet.from_json(query_path)
try:
    queries = preprocess_colabfold_msas(
        inference_query_set=queries,
        compute_settings=settings,
        use_templates=False,
    )
    (output / "query_msa.json").write_text(queries.model_dump_json(indent=4))
except MsaTransportFailure as exc:
    (output / "failure.json").write_text(json.dumps(exc.detail, ensure_ascii=False, indent=2))
    print(json.dumps(exc.detail, ensure_ascii=False), flush=True)
    raise SystemExit(1)
finally:
    (output / "network.json").write_text(json.dumps(transport.report(), ensure_ascii=False, indent=2))
    settings.cleanup_workspace()
print(json.dumps({"query_file": str(output / "query_msa.json")}), flush=True)
