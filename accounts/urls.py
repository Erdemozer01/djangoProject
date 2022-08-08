from django.urls import path, include
from .views import UsersView, MessageDetail, MessageDeleteView, UserEditView, UserAddView, UserDeleteView, \
    ProfileView, UserUpdateView, ProfileUpdateView, PasswordChance, blog_dashboard, ProfileDetailView, AddBottomView, \
    delete_bottom, AddAboutView, delete_about, AddTermsView, delete_terms, AddTitlesView, delete_titles, PostsDashBoardView
from accounts.views import UserRegister

urlpatterns = [
    path('', include('django.contrib.auth.urls')),
    path('register/', UserRegister.as_view(), name="register"),
    path('dashboard/users/', UsersView.as_view(), name="users_dashboard"),
    path('dashboard/posts/', PostsDashBoardView.as_view(), name="posts_dashboard"),
    path('blog/models/', blog_dashboard, name="blog_dashboard"),
    path('bottom/models/', AddBottomView.as_view(), name="add_bottom_dashboard"),
    path('bottom/delete/', delete_bottom, name="delete_bottom_dashboard"),
    path('about/delete/', delete_about, name="delete_about_dashboard"),
    path('about/models/', AddAboutView.as_view(), name="add_about_dashboard"),
    path('terms/models/', AddTermsView.as_view(), name="add_terms_dashboard"),
    path('title/models/', AddTitlesView.as_view(), name="add_title_dashboard"),
    path('terms/delete/', delete_terms, name="delete_terms_dashboard"),
    path('title/delete/', delete_titles, name="delete_title_dashboard"),
    path('profile/detail/<int:pk>/<slug:user>/', ProfileDetailView.as_view(), name="profile_detay"),
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
