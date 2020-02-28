
from waitress import serve
import apps
serve(apps.app, host='0.0.0.0', port=8080)
