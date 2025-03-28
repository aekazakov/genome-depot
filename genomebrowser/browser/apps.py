from django.apps import AppConfig

class BrowserConfig(AppConfig):
    name = 'browser'

    def ready(self):
        # everytime server restarts
        import browser.signals # noqa: F401

