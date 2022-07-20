from django import template
from blog.models import Cover, Bottom

register = template.Library()


@register.simple_tag()
def cover():
    return Cover.objects.all()


@register.simple_tag()
def bottom():
    return Cover.objects.all()