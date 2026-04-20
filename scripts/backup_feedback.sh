#!/bin/sh
set -eu

: "${BACKUP_REMOTE:?BACKUP_REMOTE must be set (e.g. gdrive:bottle-feedback-backups)}"
: "${FEEDBACK_DB:=/feedback/feedback.db}"

if [ ! -f "$FEEDBACK_DB" ]; then
    echo "No DB at $FEEDBACK_DB; nothing to back up."
    exit 0
fi

STAMP=$(date -u +%Y%m%dT%H%M%SZ)
SNAPSHOT="/tmp/feedback-${STAMP}.db"

sqlite3 "$FEEDBACK_DB" ".backup '$SNAPSHOT'"
rclone copyto "$SNAPSHOT" "${BACKUP_REMOTE}/feedback-${STAMP}.db"
rm -f "$SNAPSHOT"

rclone delete --min-age 30d "$BACKUP_REMOTE" || true

echo "Backup complete: ${BACKUP_REMOTE}/feedback-${STAMP}.db"
