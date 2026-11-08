# Marine Agent Deployment

This directory describes the Linux runtime used for the marine scientific
agent API.

- Runtime: `python:3.12-slim` in Docker
- Internal listener: `100.64.0.1:8098` on `lebb-node-1`
- Public paths: `/marine-agent/` and `/doho/research-api/` on `bce.dohoo.chat`
- State: `/opt/marine-agent/app/results/marine-agent-api/state.sqlite3`
- Service: `marine-agent.service`
- OpenFold3 execution: the Linux API is the control plane; Apple-MPS
  inference runs on the registered Mac worker at `100.64.0.22:8099` over
  Tailscale. The worker returns a bounded artifact bundle, which the API
  stores under `results/openfold3-remote/` for the viewer and report APIs.
- AutoDock Vina execution: the Linux API is also only the control plane;
  validated PDBQT jobs run on the same Mac at `100.64.0.22:8100`. The Vina
  worker uses a bounded parallel pool for independent ligands and returns
  artifacts under the run ID. The API image intentionally contains no Vina
  binary.

The image build context must contain `app/` (the repository runtime source)
and `wheels/` (Linux wheels downloaded for the target Python version). Model
credentials belong in `/etc/marine-agent.env` with mode `0600`; they must not
be placed in the frontend bundle or this repository.

The Nginx locations use a separate Basic Auth file while the API is being
integrated with the existing DOHO login. Anonymous requests therefore receive
`401`, and existing site paths remain unchanged.

The worker must be running before submitting a structure task. Start it on the
Mac with:

```sh
exec .env/bin/python3 openfold3/worker.py --host 100.64.0.22 --port 8099 \
  --token-file ~/.config/marine-agent/openfold-worker-token
```

The token file is local secret material and is not part of this repository.

The Vina worker is installed as `chat.marine.vina-worker` and uses the separate
`~/.config/marine-agent/vina-worker-token` file. It runs with a CPU budget
matching the Mac and does not expose a shell or arbitrary path access.

## Guarded release tools

`publish_backend.py` runs on the Linux API host and accepts a release directory,
explicit old/new images and a passed verification receipt. It checks exact file
hashes, active scientific jobs and the active image, backs up SQLite and host
files, and rolls back on health or artifact verification failure. Use
`--dry-run` for a read-only preflight. A stable high-priority systemd override
avoids accumulating version-specific override ordering bugs.

`publish_frontend.py` runs on the web host. It accepts a release directory and
the receipt for the currently reviewed frontend, preserves the homepage and
`current` symlink, and atomically replaces only declared marine-agent files.
Use `--dry-run` before publication.

Both tools require explicit operator invocation. They do not fetch code, read
private environment files into logs, commit source, or start scientific work.
