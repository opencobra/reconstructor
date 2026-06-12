from django.db import models
from django.core.validators import MinValueValidator, MaxValueValidator

class User(models.Model):
    name = models.CharField(max_length=100, blank=True, null=True)
    password = models.CharField(max_length=128, blank=True, null=True)
    saved_reactions = models.ManyToManyField(
        'Reaction', blank=True, related_name='saved_by_users')
    cred_add_to_vmh = models.BooleanField(default=False)
    cred_add_to_rhea = models.BooleanField(default=False)
    orchid_id = models.CharField(max_length=255, blank=True, null=True)
    email = models.EmailField(blank=True, null=True)
    full_name = models.CharField(max_length=255, blank=True, null=True)
    templates = models.ManyToManyField(
        'ReactionTemplate',
        blank=True,
        related_name='created_by_users')

    def check_password(self, raw_password):
        return self.password == raw_password


class Reaction(models.Model):
    substrates = models.TextField(
        help_text='Comma-separated list of substrates.')
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
    Organs = models.TextField(blank=True, null=True)
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
    stereo_counts = models.TextField(blank=True, null=True)
    stereo_locations_list = models.TextField(blank=True, null=True)
    flags = models.ManyToManyField(
        'Flag', blank=True, related_name='flagged_reactions')
    reaction_signature = models.TextField(blank=True, null=True, help_text="Unique identifier for reaction matching")
    description = models.TextField(blank=True, null=True) 

    metabolite_smiles = models.TextField(blank=True, null=True) 
    metabolite_inchis = models.TextField(blank=True, null=True)
    metabolite_inchi_keys = models.TextField(blank=True, null=True)
    metabolite_mol_weights = models.TextField(blank=True, null=True)

class ReactionsAddedVMH(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    user_name = models.CharField(max_length=255)
    user_full_name = models.CharField(max_length=255, blank=True, default='')
    user_email = models.EmailField(blank=True, default='')
    local_reaction = models.ForeignKey(
        Reaction,
        on_delete=models.SET_NULL,
        blank=True,
        null=True,
        related_name='vmh_addition_records')
    reaction_id = models.CharField(max_length=255)
    reaction_formula = models.TextField()
    reaction_abbr = models.CharField(max_length=255)
    reaction_name = models.TextField(blank=True, default='')
    direction = models.CharField(max_length=64, blank=True, default='')
    subsystem = models.TextField(blank=True, default='')
    confidence_score = models.CharField(max_length=10, blank=True, default='')
    vmh_database = models.CharField(max_length=255, blank=True, default='')
    added_metabolites = models.JSONField(blank=True, default=list)
    substrate_snapshot = models.JSONField(blank=True, default=list)
    product_snapshot = models.JSONField(blank=True, default=list)
    reference_snapshot = models.JSONField(blank=True, default=list)
    external_link_snapshot = models.JSONField(blank=True, default=list)
    gene_info_snapshot = models.JSONField(blank=True, default=list)
    comment_snapshot = models.JSONField(blank=True, default=list)
    request_snapshot = models.JSONField(blank=True, default=dict)
    matlab_result = models.JSONField(blank=True, default=dict)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.reaction_abbr} by {self.user.name}"


class MetabolitesAddedVMH(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    user_name = models.CharField(max_length=255)
    user_full_name = models.CharField(max_length=255, blank=True, default='')
    user_email = models.EmailField(blank=True, default='')
    reaction_entry = models.ForeignKey(
        ReactionsAddedVMH,
        on_delete=models.SET_NULL,
        blank=True,
        null=True,
        related_name='metabolite_records')
    local_reaction = models.ForeignKey(
        Reaction,
        on_delete=models.SET_NULL,
        blank=True,
        null=True,
        related_name='vmh_metabolite_addition_records')
    reaction_abbr = models.CharField(max_length=255, blank=True, default='')
    reaction_formula = models.TextField(blank=True, default='')
    metabolite_id = models.CharField(max_length=255)
    metabolite_formula = models.TextField(blank=True, default='')
    metabolite_abbr = models.CharField(max_length=255)
    metabolite_name = models.TextField(blank=True, default='')
    side = models.CharField(max_length=32, blank=True, default='')
    compartment = models.CharField(max_length=32, blank=True, default='')
    stoichiometry = models.CharField(max_length=64, blank=True, default='')
    source_type = models.CharField(max_length=64, blank=True, default='')
    source_identifier = models.TextField(blank=True, default='')
    inchi_key = models.CharField(max_length=255, blank=True, default='')
    vmh_database = models.CharField(max_length=255, blank=True, default='')
    submission_snapshot = models.JSONField(blank=True, default=dict)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.metabolite_abbr} by {self.user.name}"


class Subsystem(models.Model):
    name = models.CharField(max_length=255, unique=True)

    def __str__(self):
        return self.name


class CreatedReaction(models.Model):
    user = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        related_name='created_reactions')
    reaction = models.ForeignKey(
        Reaction,
        on_delete=models.CASCADE,
        related_name='created_by_users')
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.reaction.short_name} by {self.user.name}"


class Flag(models.Model):
    # Updated field name and made it optional
    name_flag = models.CharField(max_length=255, blank=True, null=True)
    color = models.CharField(max_length=7)  # Hex color code like '#FFFF00'
    user = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        related_name='flags')
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return f"{self.name_flag or 'Unnamed Flag'} ({self.color}) by {self.user.name}"


class ReactionTemplate(models.Model):
    """
    Stores user-created (or default) reaction templates
    that mirror key fields in the Reaction model.
    """
    name = models.CharField(max_length=255)
    user = models.ForeignKey(
        'User',  # or your custom user model
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        help_text="User who created this template. Null if it is a default template."
    )
    is_default = models.BooleanField(
        default=False,
        help_text="True for built-in templates; these shouldn't be edited or deleted.")

    substrates = models.TextField(
        blank=True,
        null=True,
        help_text='Comma-separated list of substrates.'
    )
    products = models.TextField(
        blank=True,
        null=True,
        help_text='Comma-separated list of products.'
    )
    direction = models.TextField(
        blank=True,
        null=True,
        help_text='Reaction direction (e.g. "forward", "bidirectional")'
    )

    substrates_types = models.TextField(blank=True, null=True)
    products_types = models.TextField(blank=True, null=True)

    subsystem = models.TextField(blank=True, null=True)
    subs_comps = models.TextField(blank=True, null=True)
    prods_comps = models.TextField(blank=True, null=True)
    subs_sch = models.TextField(blank=True, null=True)
    prods_sch = models.TextField(blank=True, null=True)
    Organs = models.TextField(
        blank=True,
        null=True
    )
    description = models.TextField(blank=True, null=True)
    def __str__(self):
        return f"{self.name} (Default: {self.is_default})"

class SavedMetabolite(models.Model):
    owner = models.ForeignKey(
        User, 
        on_delete=models.CASCADE,
        related_name='owned_metabolites'
    )
    shared_with = models.ManyToManyField(
        User,
        related_name='shared_metabolites',
        blank=True
    )
    name = models.CharField(max_length=255) # Common name
    inchi_key = models.CharField(max_length=27)  # InChIKey is 27 characters
    inchi = models.TextField(blank=True, null=True) # IUPAC International Chemical Identifier
    smiles = models.TextField(blank=True, null=True) # Simplified molecular-input line-entry system
    mol_w = models.FloatField(blank=True, null=True) # Molecular weight
    mol_file = models.TextField()  # Store the MDL molfile
    mol_formula = models.TextField(blank=True, null=True) # Molecular formula
    vmh_abbr = models.CharField(max_length=255, blank=True, null=True) # VMH abbreviation
    source_type = models.CharField(max_length=50, choices=[ 
        ('draw', 'Draw'),
        ('mol_file', 'MDL Mol File'),
        ('chebi', 'ChEBI'),
        ('swisslipids', 'SwissLipids'),
        ('pubchem', 'PubChem'),
    ])
    original_identifier = models.TextField(blank=True, null=True)  # Changed from CharField
    date_created = models.DateTimeField(auto_now_add=True) # Date the metabolite was saved

        # External Links Section
    keggId = models.CharField(max_length=50, null=True, blank=True)
    pubChemId = models.CharField(max_length=50, null=True, blank=True)
    cheBlId = models.CharField(max_length=50, null=True, blank=True)
    hmdb = models.CharField(max_length=50, null=True, blank=True)
    chembl = models.CharField(max_length=50, null=True, blank=True)
    biggId = models.CharField(max_length=191, null=True, blank=True)
    lmId = models.CharField(max_length=150, null=True, blank=True)
    ehmnId = models.CharField(max_length=50, null=True, blank=True)
    hepatonetId = models.CharField(max_length=50, null=True, blank=True)
    pdmapName = models.CharField(max_length=150, null=True, blank=True)
    biocyc = models.CharField(max_length=150, null=True, blank=True)
    chemspider = models.CharField(max_length=150, null=True, blank=True)
    drugbank = models.CharField(max_length=150, null=True, blank=True)
    food_db = models.CharField(max_length=150, null=True, blank=True)
    wikipedia = models.CharField(max_length=150, null=True, blank=True)
    metanetx = models.CharField(max_length=50, null=True, blank=True)
    seed = models.CharField(max_length=50, null=True, blank=True)
    knapsack = models.CharField(max_length=150, null=True, blank=True)
    metlin = models.CharField(max_length=150, null=True, blank=True)
    casRegistry = models.CharField(max_length=150, null=True, blank=True)
    iupac = models.TextField(null=True, blank=True)
    epa_id = models.CharField(max_length=150, null=True, blank=True)
    echa_id = models.CharField(max_length=150, null=True, blank=True)
    fda_id = models.CharField(max_length=150, null=True, blank=True)
    iuphar_id = models.CharField(max_length=150, null=True, blank=True)
    mesh_id = models.CharField(max_length=150, null=True, blank=True)
    chodb_id = models.CharField(max_length=150, null=True, blank=True)

    def __str__(self):
        return f"{self.name} ({self.inchi_key})"
    
class Workspace(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE, related_name='workspace')
    reactions = models.ManyToManyField(Reaction, related_name='in_workspaces')
    updated_at = models.DateTimeField(auto_now=True)

    def __str__(self):
        return f"{self.user.name}'s Workspace"


class ReactionGroup(models.Model):
    """
    A named, per-user, ordered collection of reactions shown under one header
    tab in the multi-reaction workspace (feature 001-multi-reaction-tabbed).

    Distinct from `Workspace` (one-per-user VMH-staging bag): a user may own many
    named groups. Membership ordering is carried by `ReactionGroupMembership`.
    Deleting a group does NOT delete the underlying `Reaction` rows it referenced.
    """
    owner = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        related_name='reaction_groups')
    name = models.CharField(max_length=255)
    # Exactly one group per user should be active (last-opened tab).
    is_active = models.BooleanField(default=False)
    reactions = models.ManyToManyField(
        Reaction,
        through='ReactionGroupMembership',
        related_name='reaction_groups')
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        ordering = ['created_at']

    def __str__(self):
        return f"{self.name} (group of {self.owner.name})"


class ReactionGroupMembership(models.Model):
    """Ordered through-model linking a `Reaction` into a `ReactionGroup`."""
    group = models.ForeignKey(
        ReactionGroup,
        on_delete=models.CASCADE,
        related_name='memberships')
    reaction = models.ForeignKey(
        Reaction,
        on_delete=models.CASCADE,
        related_name='group_memberships')
    position = models.PositiveIntegerField(default=0)
    added_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ['position', 'added_at']
        # A reaction may appear at most once per group (duplicate-guard, FR-012).
        unique_together = ('group', 'reaction')

    def __str__(self):
        return f"{self.reaction_id} in {self.group_id} @ {self.position}"
   
# Gene table, populated by HGNC to show suggested gene name list in Gene Info interface 
class Gene(models.Model):
    symbol = models.CharField(max_length=50, db_index=True)
    name = models.TextField(blank=True)
    hgnc_id = models.CharField(max_length=20, unique=True)
    entrez_id = models.CharField(max_length=20, blank=True, null=True, db_index=True)
    aliases = models.TextField(blank=True)  # store comma-separated aliases

    def __str__(self):
        return self.symbol
