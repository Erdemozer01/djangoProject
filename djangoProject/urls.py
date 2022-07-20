from django.contrib import admin
from django.urls import path, include
from django.conf import settings
from django.conf.urls.static import static
from django.contrib.auth.views import LogoutView

urlpatterns = [
                  path('admin/', admin.site.urls),
                  path('ckeditor/', include('ckeditor_uploader.urls')),
                  path('accounts/', include('accounts.urls')),
                  path('', include('blog.urls')),
                  path('post/', include('post.urls')),
                  path('logout/', LogoutView.as_view(), name="logout"),
              ] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
