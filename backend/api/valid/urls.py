from django.urls import path

from api.valid.views import Valid

urlpatterns = [path("", Valid.as_view())]
