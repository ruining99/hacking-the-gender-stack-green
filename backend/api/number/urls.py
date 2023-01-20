from django.urls import path

from api.number.views import Number

urlpatterns = [path("", Number.as_view())]
