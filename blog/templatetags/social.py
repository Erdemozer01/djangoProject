from django import template
from blog.models import Social

register = template.Library()

@register.simple_tag()
def social():
    return Social.objects.all()

