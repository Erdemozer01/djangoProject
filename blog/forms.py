from django import forms
from blog.models import Inbox


class ContactForm(forms.ModelForm):
    terms = forms.BooleanField(required=True)

    class Meta:
        model = Inbox
        fields = '__all__'
        widgets = {
            'name': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Ad Soyad'}),
            'email': forms.EmailInput(attrs={'class': 'form-control', 'placeholder': 'Email'}),
            'topic': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Konu'}),
            'phone': forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Telefon'}),
            'message': forms.Textarea(attrs={'class': 'form-control', 'placeholder': 'Mesajınız'}),
        }
