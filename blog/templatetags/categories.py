from django import template
from post.models import Category

register = template.Library()


@register.simple_tag()
def categories():
    return Category.objects.all()

