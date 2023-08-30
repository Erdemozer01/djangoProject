from django.contrib import admin
from django.urls import path, include
from django.conf import settings
from django.conf.urls.static import static
from django.contrib.auth.views import LogoutView
from django.contrib.staticfiles.urls import staticfiles_urlpatterns

urlpatterns = [
                  path('admin/', admin.site.urls),
                  path('ckeditor/', include('ckeditor_uploader.urls')),
                  path('accounts/', include('accounts.urls')),
                  path('', include('blog.urls')),
                  path('post/', include('post.urls')),
                  path('laboratory/bioinformatic/', include('bioinformatic.urls')),
                  path('logout/', LogoutView.as_view(), name="logout"),
                  path('django_plotly_dash/', include('django_plotly_dash.urls', namespace="plotly")),
              ] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
