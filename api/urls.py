from django.urls import path
from . import views

app_name = 'main'

urlpatterns = [
    path('v1', views.v1, name='v1'),
    path('v2', views.v2, name='v2'),
    path('v3', views.v3, name='v3')
]