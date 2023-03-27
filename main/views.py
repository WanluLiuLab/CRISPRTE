from django.shortcuts import render, get_object_or_404
from django.http import HttpResponse
from django.contrib.auth.decorators import login_required
from django.template import loader
from django.conf import settings
from django import forms
from django.views.decorators.csrf import csrf_exempt
from django.http import JsonResponse
from django.views.static import serve

import os
import json

def index(request):
    template = loader.get_template('index.html')
    return HttpResponse(template.render({"version": "0.0.1 beta"}, request))

def docs(request):
    template = loader.get_template('docs.html')
    return HttpResponse(template.render({"version": "0.0.1 beta"}, request))