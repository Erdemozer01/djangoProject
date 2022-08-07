from django.urls import path
from . import views

app_name = "blog"

urlpatterns = [
    path('', views.HomeView.as_view(), name='home'),
    path('terms/', views.terms, name='terms'),
    path('hakkımızda/', views.about, name='about'),
    path('iletisim/', views.ContactView.as_view(), name='contact'),
]
