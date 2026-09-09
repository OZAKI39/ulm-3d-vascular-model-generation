/* Reuse the existing offline Chrome/CDP lifecycle, interactions and cleanup. */
'use strict';
process.env.MIRHEO_REVIEW_KIND='equilibration';
require('./check_sdpd_diagnostics_browser.cjs');
