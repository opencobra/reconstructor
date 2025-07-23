from django.apps import AppConfig
from reactions.utils.MatlabSessionManager import MatlabSessionManager


class ReactionsConfig(AppConfig):
    default_auto_field = "django.db.models.BigAutoField"
    name = "reactions"
    
   
