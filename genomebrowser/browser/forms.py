from django import forms
from django.forms.widgets import TextInput
from browser.models import Tag
from browser.models import Config


class TsvImportForm(forms.Form):
    tsv_file = forms.FileField()


class ExcelImportForm(forms.Form):
    xlsx_file = forms.FileField()


class GenomeImportForm(forms.Form):
    tsv_file = forms.FileField()
    zip_file = forms.FileField(required=False,
                               label='Zip archive with genomes in ' + 
                               'GBFF format (optional)'
                               )
    download_email = forms.EmailField(max_length=200,
                                      required=False,
                                      label='Email for NCBI genome downloads (optional)'
                                      )


class TagModelForm(forms.ModelForm):
    class Meta:
        model = Tag
        fields = "__all__"
        widgets = {
            "color": TextInput(attrs={"type": "color"}),
            "textcolor": TextInput(attrs={"type": "color"}),
        }


class AddTagForm(forms.Form):
    tag = forms.ModelChoiceField(queryset=Tag.objects.all(), empty_label='Not selected')


class RemoveTagForm(forms.Form):
    tag = forms.ModelChoiceField(queryset=Tag.objects.all(), empty_label='Not selected')


class ChooseAnnotationToolForm(forms.Form):
    choices = []
    for item in Config.objects.values_list('param', 'value'):
        if item[0].startswith('plugins.') and item[0].endswith('.display_name'):
            choices.append((item[0].split('.')[1], item[1]))
    tools = forms.MultipleChoiceField(
        choices = choices,
        widget=forms.CheckboxSelectMultiple,
    )
    
