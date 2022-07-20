from django import template
from accounts.models import Profile

register = template.Library()


@register.simple_tag()
def profile():
    return Profile.objects.all()
