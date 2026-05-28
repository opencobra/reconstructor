from django import forms
from .models import Reaction


class ReactionForm(forms.ModelForm):
    class Meta:
        model = Reaction
        fields = ['substrates', 'products']

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['substrates'].required = False
        self.fields['products'].required = False

    def clean_substrates(self):
        data = self.cleaned_data['substrates']
        # You can add more data cleaning or validation here if needed
        return data

    def clean_products(self):
        data = self.cleaned_data['products']
        # You can add more data cleaning or validation here if needed
        return data
