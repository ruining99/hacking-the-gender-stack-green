from django.urls import path

from api.lipinski.views import Lipinski

urlpatterns = [path("", Lipinski.as_view())]
