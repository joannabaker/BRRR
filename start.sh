export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/library}"
mkdir -p "$R_LIBS_USER"

#!/usr/bin/env bash
set -euo pipefail

# Default port 3838 unless you pass one: ./start.sh 4000
PORT="${1:-3838}"

R -q -e "shiny::runApp('app.R', host='0.0.0.0', port=${PORT}, launch.browser=FALSE)"
