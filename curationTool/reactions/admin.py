from django.contrib import admin
from .models import User, Reaction, ReactionsAddedVMH, MetabolitesAddedVMH, Subsystem, CreatedReaction

class CreatedReactionInline(admin.TabularInline):
    model = CreatedReaction
    extra = 0  # No extra empty forms displayed
    fields = ('reaction', 'created_at')  # Fields to display
    readonly_fields = ('created_at',)  # Make created_at read-only

class UserAdmin(admin.ModelAdmin):
    list_display = ('id', 'name', 'cred_add_to_vmh', 'cred_add_to_rhea')  # You can add more fields here
    search_fields = ('id', 'name')  # Fields to search by in the admin site
    filter_horizontal = ('saved_reactions',)  # For many-to-many fields
    inlines = [CreatedReactionInline]  # Add inline admin for CreatedReaction

class ReactionAdmin(admin.ModelAdmin):
    list_display = ('id', 'short_name', 'references', 'ext_links', 'gene_info', 'comments')  # You can add more fields here
    search_fields = ('id', 'short_name')  # Fields to search by in the admin site

class CreatedReactionAdmin(admin.ModelAdmin):
    list_display = ('id', 'user', 'reaction', 'created_at')  # Customize list display
    search_fields = ('user__name', 'reaction__short_name')  # Customize search fields

admin.site.register(User, UserAdmin)
admin.site.register(Reaction, ReactionAdmin)
admin.site.register(ReactionsAddedVMH)
admin.site.register(MetabolitesAddedVMH)
admin.site.register(Subsystem)
admin.site.register(CreatedReaction, CreatedReactionAdmin)
