from django.urls import path
from . import views

app_name = 'main'

urlpatterns = [
    path('', views.index, name='crisprte'),
    path('doc', views.docs, name='doc')
]