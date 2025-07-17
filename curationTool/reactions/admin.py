from django.contrib import admin
from .models import User, Reaction, ReactionsAddedVMH, MetabolitesAddedVMH, Subsystem, CreatedReaction, Flag, ReactionTemplate, SavedMetabolite, Workspace



class CreatedReactionInline(admin.TabularInline):
    model = CreatedReaction
    extra = 0  # No extra empty forms displayed
    fields = ('reaction', 'created_at')  # Fields to display
    readonly_fields = ('created_at',)  # Make created_at read-only


class FlagInline(admin.TabularInline):
    model = Flag
    extra = 1  # Allows adding one empty form by default
    fields = ('name_flag', 'color', 'created_at')  # Fields to display
    readonly_fields = ('created_at',)  # Make created_at read-only


class UserAdmin(admin.ModelAdmin):
    list_display = (
        'id',
        'name',
        'cred_add_to_vmh',
        'cred_add_to_rhea')  # You can add more fields here
    search_fields = ('id', 'name')  # Fields to search by in the admin site
    filter_horizontal = ('saved_reactions',)  # For many-to-many fields
    inlines = [CreatedReactionInline]  # Add inline admin for CreatedReaction


class ReactionAdmin(admin.ModelAdmin):
    list_display = (
        'id',
        'short_name',
        'references',
        'ext_links',
        'gene_info',
        'comments')  # You can add more fields here
    # Fields to search by in the admin site
    search_fields = ('id', 'short_name')


class CreatedReactionAdmin(admin.ModelAdmin):
    list_display = (
        'id',
        'user',
        'reaction',
        'created_at')  # Customize list display
    # Customize search fields
    search_fields = ('user__name', 'reaction__short_name')


class FlagAdmin(admin.ModelAdmin):
    # Customize list display for flags
    list_display = ('id', 'name_flag', 'color', 'user', 'created_at')
    search_fields = (
        'name_flag',
        'user__name',
        'color')  # Customize search fields
    list_filter = ('color', 'user')  # Add filters for color and user


class ReactionTemplateAdmin(admin.ModelAdmin):
    """
    Admin configuration for ReactionTemplate model.
    """
    list_display = (
        'id',
        'name',
        'is_default',
        'user',
        'direction',
        'subsystem')  # Customize list display
    search_fields = ('name', 'user__name', 'subsystem',
                     'direction')  # Fields to search by
    # Add filters for default templates, users, and subsystems
    list_filter = ('is_default', 'user', 'subsystem')
    readonly_fields = ('is_default',)  # Make is_default read-only
    fieldsets = (
        (None, {
            'fields': ('name', 'user', 'is_default')
        }),
        ('Reaction Details', {
            'fields': ('direction', 'subsystem', 'Organs')
        }),
        ('Substrates', {
            'fields': ('substrates', 'substrates_types', 'subs_comps', 'subs_sch')
        }),
        ('Products', {
            'fields': ('products', 'products_types', 'prods_comps', 'prods_sch')
        }),
    )

class SavedMetaboliteAdmin(admin.ModelAdmin):
    list_display = (
        'id',
        'owner',
        'name',
        'inchi_key',
        'vmh_abbr',
        'source_type',
        'original_identifier',
        'date_created')  # Customize list display
    search_fields = ('name', 'inchi_key', 'vmh_abbr',
                     'source_type', 'original_identifier')  # Fields to search by
    list_filter = ('source_type',)  # Add filter for source_type
    readonly_fields = ('date_created',)  # Make date_created read-only
admin.site.register(User, UserAdmin)
admin.site.register(Reaction, ReactionAdmin)
admin.site.register(ReactionsAddedVMH)
admin.site.register(MetabolitesAddedVMH)
admin.site.register(Subsystem)
admin.site.register(CreatedReaction, CreatedReactionAdmin)
admin.site.register(Flag, FlagAdmin)
admin.site.register(ReactionTemplate, ReactionTemplateAdmin)
admin.site.register(SavedMetabolite, SavedMetaboliteAdmin)
admin.site.register(Workspace)  