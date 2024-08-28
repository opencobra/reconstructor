from django.db import models

class User(models.Model):
    name = models.CharField(max_length=100, blank=True, null=True)
    password = models.CharField(max_length=128, blank=True, null=True)
    saved_reactions = models.ManyToManyField('Reaction', blank=True, related_name='saved_by_users')
    cred_add_to_vmh = models.BooleanField(default=False)
    cred_add_to_rhea = models.BooleanField(default=False)
    orchid_id = models.CharField(max_length=255, blank=True, null=True)
    email = models.EmailField(blank=True, null=True)
    full_name = models.CharField(max_length=255, blank=True, null=True)

    def check_password(self, raw_password):
        return raw_password == self.password


class Reaction(models.Model):
    substrates = models.TextField(help_text='Comma-separated list of substrates.')
    products = models.TextField(help_text='Comma-separated list of products.')
    short_name = models.TextField(blank=True, null=True)
    direction = models.TextField(blank=True, null=True)
    substrates_types = models.TextField(blank=True, null=True)
    products_types = models.TextField(blank=True, null=True)
    substrates_names = models.TextField(blank=True, null=True)
    products_names = models.TextField(blank=True, null=True)
    visualization = models.TextField(blank=True, null=True)
    molc_formula = models.TextField(blank=True, null=True)
    balanced_count = models.TextField(blank=True, null=True)
    balanced_charge = models.TextField(blank=True, null=True)
    subsystem = models.TextField(blank=True, null=True)
    subs_comps = models.TextField(blank=True, null=True)
    prods_comps = models.TextField(blank=True, null=True)
    subs_sch = models.TextField(blank=True, null=True)
    prods_sch = models.TextField(blank=True, null=True)
    subs_atoms = models.TextField(blank=True, null=True)
    prods_atoms = models.TextField(blank=True, null=True)
    subs_charge = models.TextField(blank=True, null=True)
    prods_charge = models.TextField(blank=True, null=True)
    symb_to_name = models.TextField(blank=True, null=True)
    metabolite_names = models.TextField(blank=True, null=True)
    metabolite_formulas = models.TextField(blank=True, null=True)
    metabolite_charges = models.TextField(blank=True, null=True)
    metabolite_mol_file_strings = models.TextField(blank=True, null=True)
    subs_found = models.TextField(blank=True, null=True)
    subs_miriams = models.TextField(blank=True, null=True)
    prod_found = models.TextField(blank=True, null=True)
    prod_miriams = models.TextField(blank=True, null=True)
    Organs = models.TextField(blank=True, null=True)    # New field
    vmh_found = models.BooleanField(default=False)
    vmh_found_similar = models.BooleanField(default=False)
    vmh_url = models.TextField(blank=True, null=True)
    vmh_formula = models.TextField(blank=True, null=True)
    references = models.JSONField(blank=True, null=True)
    ext_links = models.JSONField(blank=True, null=True)
    gene_info = models.JSONField(blank=True, null=True)
    comments = models.JSONField(blank=True, null=True)
    confidence_score = models.CharField(blank=True, max_length=10, null=True)
    rxn_formula = models.TextField(blank=True, null=True)

    flags = models.ManyToManyField('Flag', blank=True, related_name='flagged_reactions')


class ReactionsAddedVMH(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    user_name = models.CharField(max_length=255)
    reaction_id = models.CharField(max_length=255)
    reaction_formula = models.TextField()
    reaction_abbr = models.CharField(max_length=255)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.reaction_abbr} by {self.user.name}"


class MetabolitesAddedVMH(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    user_name = models.CharField(max_length=255)
    metabolite_id = models.CharField(max_length=255)
    metabolite_formula = models.TextField()
    metabolite_abbr = models.CharField(max_length=255)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.metabolite_abbr} by {self.user.name}"


class Subsystem(models.Model):
    name = models.CharField(max_length=255, unique=True)

    def __str__(self):
        return self.name


class CreatedReaction(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name='created_reactions')
    reaction = models.ForeignKey(Reaction, on_delete=models.CASCADE, related_name='created_by_users')
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.reaction.short_name} by {self.user.name}"

class Flag(models.Model):
    name_flag = models.CharField(max_length=255, blank=True, null=True)  # Updated field name and made it optional
    color = models.CharField(max_length=7)  # Hex color code like '#FFFF00'
    user = models.ForeignKey(User, on_delete=models.CASCADE, related_name='flags')
    created_at = models.DateTimeField(auto_now_add=True)
    
def __str__(self):
    return f"{self.name_flag or 'Unnamed Flag'} ({self.color}) by {self.user.name}"