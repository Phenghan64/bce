FROM python:3.12-slim

WORKDIR /app

# Scientific execution is delegated to the authenticated Mac worker.  The
# control-plane image intentionally does not contain a Vina binary; a worker
# outage is reported as a structured task failure instead of silently falling
# back to the API node.

# The deployment context supplies the application source and an offline wheel
# directory. Keeping dependency installation in the image makes the runtime
# independent from the host Python installation.
COPY app/ /app/
COPY wheels/ /opt/marine-agent-wheels/

RUN python -m pip install --no-cache-dir --no-index \
      --find-links=/opt/marine-agent-wheels \
      -r /app/marine_agent/requirements.txt \
    && rm -rf /opt/marine-agent-wheels

ENV PYTHONPATH=/app
ENV PYTHONUNBUFFERED=1

EXPOSE 8098

CMD ["uvicorn", "marine_agent.api:app", "--host", "0.0.0.0", "--port", "8098"]
