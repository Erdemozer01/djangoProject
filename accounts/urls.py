"""myblog URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.urls import path, include
from .views import DashBoardView, MessageDetail, MessageDeleteView, UserEditView, UserAddView, UserDeleteView, \
    ProfileView, edit, UserUpdateView, ProfileUpdateView, PasswordChance
from accounts.views import UserRegister

urlpatterns = [
    path('register/', UserRegister.as_view(), name="register"),
    path('dashboard/', DashBoardView.as_view(), name="dashboard"),
    path('edit/<int:pk>/<slug:username>/', UserEditView.as_view(), name="useredit"),
    path('profile/<int:pk>/<slug:username>/', ProfileView.as_view(), name="profile"),
    path('user_edit/<int:pk>/<slug:username>/', UserUpdateView.as_view(), name="user_update"),
    path('profile_edit/<pk>/<username>/', ProfileUpdateView.as_view(), name="profile_update"),
    path('add_user/', UserAddView.as_view(), name="user_add"),
    path('add_user/<int:pk>/', UserDeleteView.as_view(), name="user_del"),
    path('messages/<int:pk>/<slug:name>/', MessageDetail.as_view(), name="messages"),
    path('messages/<int:pk>/', MessageDeleteView.as_view(), name="message_delete"),
    path('password/<pk>/<username>/', PasswordChance.as_view(), name="password_change"),
]
