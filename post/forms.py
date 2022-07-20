from django import forms
from .models import Posts


class PostForm(forms.ModelForm):
    class Meta:
        model = Posts
        fields = "__all__"
        labels = {
            'text': ""
        }

        widgets = {
            'cover': forms.FileInput(attrs={'class': 'form-control'}),
            'title': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Başlık'}),
            'tags': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Etiket'}),
            'author': forms.Select(attrs={'class': 'form-control'}),
            'category': forms.Select(attrs={'class': 'form-control'}),
            'publish': forms.DateInput(attrs={'class': 'form-control'}),
            'status': forms.Select(attrs={'class': 'form-control'})
        }
