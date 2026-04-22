#!/usr/bin/env bash
set -euo pipefail

required_vars=(
  AUTH_REDIRECT_URI
  AUTH_COOKIE_SECRET
  AUTH_GOOGLE_CLIENT_ID
  AUTH_GOOGLE_CLIENT_SECRET
  AUTH_GOOGLE_SERVER_METADATA_URL
)
for v in "${required_vars[@]}"; do
  if [[ -z "${!v:-}" ]]; then
    echo "entrypoint: required env var $v is not set" >&2
    exit 1
  fi
done

mkdir -p /app/.streamlit
umask 077
cat > /app/.streamlit/secrets.toml <<EOF
[auth]
redirect_uri = "${AUTH_REDIRECT_URI}"
cookie_secret = "${AUTH_COOKIE_SECRET}"

[auth.google]
client_id = "${AUTH_GOOGLE_CLIENT_ID}"
client_secret = "${AUTH_GOOGLE_CLIENT_SECRET}"
server_metadata_url = "${AUTH_GOOGLE_SERVER_METADATA_URL}"
EOF

exec streamlit run path_viewer/app.py \
  --server.address=0.0.0.0 \
  --server.port="${PORT:-8501}" \
  --theme.base="light"
