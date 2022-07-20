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
from django.urls import path
from . import views

app_name = "post"

urlpatterns = [
    path('categories/', views.CategoriesView.as_view(), name='categories'),
    path('category/<slug:slug>/', views.CategoryView.as_view(), name='category'),
    path('<slug:category>/<slug:title>/<int:pk>/', views.PostDetailView.as_view(), name='post_detail'),
    path('add/', views.AddPostView.as_view(), name='add_post'),
    path('add_category/', views.AddCategoryView.as_view(), name='add_category'),
    path('update/<slug:slug>', views.UpdatePostView.as_view(), name='update_post'),
    path('<pk>/', views.DeletePostView.as_view(), name='delete_post'),
]
