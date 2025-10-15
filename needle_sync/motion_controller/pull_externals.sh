#!/usr/bin/env bash
set -euo pipefail

# pull_externals.sh
# Fetches header-only deps for the motion_controller stack (reproducible + shallow).
# If a target dir exists but isn't a git repo, it is backed up (*.bak.<timestamp>)
# or deleted if EXTERNALS_FORCE=1 is set.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTERNALS_DIR="$ROOT_DIR/externals"
mkdir -p "$EXTERNALS_DIR"

clone_or_update() {
  local name="$1"
  local repo_url="$2"
  local dest="$EXTERNALS_DIR/$name"

  if git -C "$dest" rev-parse --git-dir >/dev/null 2>&1; then
    git -C "$dest" remote set-url origin "$repo_url" || true
    git -C "$dest" fetch --depth=1 origin
    git -C "$dest" reset --hard FETCH_HEAD
    return
  fi

  if [[ -d "$dest" ]]; then
    if [[ "${EXTERNALS_FORCE:-}" == "1" ]]; then
      rm -rf "$dest"
    else
      ts="$(date +%Y%m%d-%H%M%S)"
      echo "Directory $dest exists but is not a git repo. Backing up to ${dest}.bak.$ts"
      mv "$dest" "${dest}.bak.$ts"
    fi
  fi

  echo "Cloning $name"
  git clone --depth=1 "$repo_url" "$dest"
}

# Ruckig: header-only OTG library.
clone_or_update "ruckig" "https://github.com/pantor/ruckig.git"

# Eigen: linear algebra.
clone_or_update "eigen" "https://gitlab.com/libeigen/eigen.git"

# Kalman: header-only EKF/UKF/SR variants, Eigen-based (replacement for OpenKalman).
clone_or_update "kalman" "https://github.com/mherb/kalman.git"

# iir1: biquad filter toolkit.
clone_or_update "iir1" "https://github.com/berndporr/iir1.git"

echo "Dependency sync complete."