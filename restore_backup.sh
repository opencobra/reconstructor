#!/bin/bash

# Restoration script for database and media files
# This script restores the PostgreSQL database and media files from backups

set -e  # Exit on error

echo "=== Starting backup restoration ==="

# Configuration
DB_BACKUP="/home/saleh/Downloads/db_backup_20250811_020001.sql"
MEDIA_BACKUP="/home/saleh/Downloads/media_backup_20250811_020001.tar"
MEDIA_DIR="/home/saleh/reconstructor/curationTool/media"
POSTGRES_DATA_DIR="/home/saleh/reconstructor/docker/postgres_data"

# Check if backup files exist
if [ ! -f "$DB_BACKUP" ]; then
    echo "Error: Database backup file not found: $DB_BACKUP"
    exit 1
fi

if [ ! -f "$MEDIA_BACKUP" ]; then
    echo "Error: Media backup file not found: $MEDIA_BACKUP"
    exit 1
fi

# Step 1: Stop any running containers
echo "Stopping existing containers..."
docker compose down

# Step 2: Reset the Postgres data directory (kept outside of Docker)
echo "Preparing Postgres data directory..."
if [ -d "$POSTGRES_DATA_DIR" ] && [ "$(ls -A "$POSTGRES_DATA_DIR")" ]; then
    echo "Backing up existing Postgres data directory..."
    mv "$POSTGRES_DATA_DIR" "${POSTGRES_DATA_DIR}.old.$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$POSTGRES_DATA_DIR"

# Step 3: Start the database service
echo "Starting database service..."
docker compose up -d db

# Wait for database to be ready
echo "Waiting for database to be ready..."
sleep 10

# Step 4: Restore the database
echo "Restoring database from backup..."
docker compose exec -T db psql -U saleh -d curationtooldb1 < "$DB_BACKUP"

if [ $? -eq 0 ]; then
    echo "✓ Database restored successfully"
else
    echo "✗ Database restoration failed"
    exit 1
fi

# Step 5: Restore media files
echo "Restoring media files..."
# Backup existing media if it exists
if [ -d "$MEDIA_DIR" ]; then
    echo "Backing up existing media directory..."
    mv "$MEDIA_DIR" "${MEDIA_DIR}.old.$(date +%Y%m%d_%H%M%S)"
fi

# Create media directory
mkdir -p "$MEDIA_DIR"

# Extract media backup with memory-efficient options
echo "Extracting media files..."
# Extract to a temporary directory first, then move the contents to media
TEMP_EXTRACT_DIR="/tmp/media_restore_$$"
mkdir -p "$TEMP_EXTRACT_DIR"
tar -xf "$MEDIA_BACKUP" -C "$TEMP_EXTRACT_DIR" --no-same-owner --checkpoint=1000 --checkpoint-action=dot

# Move the extracted directories to the media folder
if [ -d "$TEMP_EXTRACT_DIR/images" ]; then
    mv "$TEMP_EXTRACT_DIR/images" "$MEDIA_DIR/"
fi
if [ -d "$TEMP_EXTRACT_DIR/mol_files" ]; then
    mv "$TEMP_EXTRACT_DIR/mol_files" "$MEDIA_DIR/"
fi
if [ -d "$TEMP_EXTRACT_DIR/rxn_files" ]; then
    mv "$TEMP_EXTRACT_DIR/rxn_files" "$MEDIA_DIR/"
fi

# Clean up temp directory
rm -rf "$TEMP_EXTRACT_DIR"

if [ $? -eq 0 ]; then
    echo "✓ Media files restored successfully"
else
    echo "✗ Media files restoration failed"
    exit 1
fi

# Step 6: Set proper permissions
echo "Setting permissions..."
chmod -R 755 "$MEDIA_DIR"

# Step 7: Start all services
echo "Starting all services..."
docker compose up -d

echo ""
echo "=== Restoration completed successfully! ==="
echo "Your application should now be running with the restored data."
echo "Access it at: http://localhost:8001"
echo ""
