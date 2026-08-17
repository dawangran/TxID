FROM dawang02/txid:0.1.0-jupyter@sha256:c93e5ff2ca5796ad1a39669510ca5ca5c3a279ac926d56aba2f1f1b461a03439 AS builder

USER root
WORKDIR /opt/txid
RUN python -m pip install "setuptools==80.9.0"
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
COPY schemas ./schemas
RUN python -m pip wheel --no-deps --no-build-isolation --wheel-dir /tmp/wheels "."

FROM dawang02/txid:0.1.0-jupyter@sha256:c93e5ff2ca5796ad1a39669510ca5ca5c3a279ac926d56aba2f1f1b461a03439

LABEL org.opencontainers.image.title="TxID with JupyterLab" \
      org.opencontainers.image.description="Reference-aware deterministic transcript identities and JupyterLab" \
      org.opencontainers.image.source="https://github.com/dawangran/TxID" \
      org.opencontainers.image.licenses="BSD-3-Clause" \
      org.opencontainers.image.version="0.1.3"

USER root
COPY --from=builder /tmp/wheels/txid-0.1.3-py3-none-any.whl /tmp/txid-0.1.3-py3-none-any.whl
RUN python -m pip install --no-deps /tmp/txid-0.1.3-py3-none-any.whl \
    && txid --version \
    && jupyter lab --version \
    && rm /tmp/txid-0.1.3-py3-none-any.whl

USER txid
WORKDIR /workspace
EXPOSE 8888

CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--no-browser", "--ServerApp.root_dir=/workspace"]
