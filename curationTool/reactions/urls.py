from django.contrib import admin
from django.urls import path
from .views import (
    about_leaderboard_views,
    user_views,
    reaction_views,
    flag_views,
    gene_views,
    metabolite_views,
    template_views,
    vmh_views,
    ai_views,
    utility_views,
    graph_views)
from django.conf import settings
from django.conf.urls.static import static

urlpatterns = [
    path('get_user/', user_views.get_user, name='get_user'),
    path('set_session_user/', user_views.set_session_user, name='set_session_user'),
    path('register_user', user_views.register_user, name='register_user'),
    path('save_reaction/', user_views.save_user_reaction, name='save_reaction'),

    path('delete_reaction/', reaction_views.delete_reaction, name='delete_reaction'),
    path('', reaction_views.input_reaction, name='input_reaction'),
    path('add_info_to_reaction/', reaction_views.add_info_to_reaction, name='add_info_to_reaction'),
    path('get_reaction_details/', reaction_views.get_reaction_details, name='get_reaction_details'),
    path('delete_reaction_info/', reaction_views.delete_reaction_info, name='delete_reaction_info'),
    path('saved_reactions/reactions/clone/', reaction_views.clone_reaction_view, name='clone_reaction'),
    path('user-reactions-vmh/', reaction_views.get_user_reactions_and_vmh, name='get_user_reactions_and_vmh'),
    path('get_reaction/<int:reaction_id>/', reaction_views.get_reaction, name='get_reaction'),
    path('update_gene_info/', reaction_views.update_gene_info, name='update_gene_info'),
    path('identical_reaction/', reaction_views.identical_reaction, name='identical_reaction'),
    path('create-reaction/', reaction_views.create_reaction, name='create_reaction'),
    path('reaction_name_exists/', reaction_views.reaction_name_exists, name='reaction_name_exists'),
    path('check-reaction', reaction_views.already_saved, name='already_saved'),
    path('saved_reactions/', reaction_views.saved_reactions, name='saved_reactions'),
    path('available_reactions', reaction_views.get_available_reactions, name='available_reactions'),
    path('edit_reaction_info/', reaction_views.edit_reaction_info, name='edit_reaction_info'),
    path('update_confidence_score/', reaction_views.update_confidence_score, name='update_confidence_score'),

    path('saved_reactions/flags/<int:user_id>/', flag_views.get_user_flags, name='get_user_flags'),
    path('flags/<int:user_id>/', flag_views.get_user_flags, name='get_user_flags'),
    path('add_flag/', flag_views.add_flag, name='add_user_flag'),
    path('saved_reactions/add_flag/', flag_views.add_flag, name='add_user_flag'),
    path('saved_reactions/remove_flag/', flag_views.remove_flag, name='remove_user_flag'),
    path('saved_reactions/save_flags_in_saved_reactions/', flag_views.save_flags_in_saved_reactions, name='save_flags_in_saved_reactions'),

    path('get_gene_info/', gene_views.get_gene_info, name='get_gene_info'),
    path("gene_details_view/", gene_views.gene_details_view, name="gene_details_view"),
    path("parse_gene_info/", gene_views.parse_gene_info, name="parse_gene_info"),
    path('gene_parsing/', gene_views.gene_parsing, name='gene_parsing'),
    path("api/gene/suggest", gene_views.suggest_genes, name="suggest_genes"),


    path('fetch_rhea_rxn', metabolite_views.fetch_rhea_rxn, name='fetch_rhea_rxn'),
    path('verify_metabolite/', metabolite_views.verify_metabolite, name='verify_metabolite'),
    path('get_saved_metabolites/', metabolite_views.get_saved_metabolites, name='get_saved_metabolites'),
    path('update_metabolite/<int:metabolite_id>/', metabolite_views.update_metabolite, name='update_metabolite'),
    path('delete_metabolite/<int:metabolite_id>/', metabolite_views.delete_metabolite, name='delete_metabolite'),
    path('share_metabolites/', metabolite_views.share_metabolites, name='share_metabolites'),

    path('list_templates/', template_views.list_templates, name='list_templates'),
    path('create_template/', template_views.create_template, name='create_template'),
    path('share_template/', template_views.share_template, name='share_template'),
    path('update_template/', template_views.update_template, name='update_template'),
    path('delete_template/', template_views.delete_template, name='delete_template'),
    path('get_rxn_template/', template_views.get_rxn_template, name='get_rxn_template'),

    path('add_to_vmh/', vmh_views.add_to_vmh, name='add_to_vmh'),
    path('get_subsystems/', vmh_views.get_subsystems, name='get_subsystems'),
    path('prepare_add_to_vmh/', vmh_views.prepare_add_to_vmh, name='prepare_add_to_vmh'),
    path('check_vmh_availability/', vmh_views.check_vmh_availability, name='check_vmh_availability'),
    path('update_subsystems/', vmh_views.update_subsystems, name='update_subsystems'),
    path('create-formula-abbr/', vmh_views.create_formula_abbr, name='create_formula_abbr'),
    path('check_reaction_vmh/', vmh_views.check_reaction_vmh, name='check_reaction_vmh'),
    path('get_from_vmh/', vmh_views.get_from_vmh, name='get_from_vmh'),
    
    path('send_to_workspace/', vmh_views.send_to_workspace, name='send_to_workspace'),
    path('remove_from_workspace/', vmh_views.remove_from_workspace, name='remove_from_workspace'),
    path('VMH_Workspace/', vmh_views.vmh_workspace, name='VMH_Workspace'),
    path('save_reaction_draft/', vmh_views.save_reaction_draft, name='save_reaction_draft'),

    path('get_ai_response/', ai_views.get_ai_response, name='get_ai_response'),

    path('stats/', about_leaderboard_views.leader_board, name='leader_board'),
    path('about/', about_leaderboard_views.about_view, name='about'),

    path('check-session/', utility_views.check_session_data, name='check_session'),
    path('delete-gene-info/', utility_views.delete_gene_info_from_session, name='delete_gene_info_from_session'),
    path('clear-session/', utility_views.clear_session, name='clear_session'),
    path('save-formula/', utility_views.save_formula, name='save_formula'),
    path('parse_formula_with_compartments/', utility_views.parse_formula_with_compartments, name='parse_formula_with_compartments'),
    path('get_pubmed_info/<str:pmid>/', utility_views.get_pubmed_info, name='get_pubmed_info'),
    path('chemdoodle_sketcher/', utility_views.chemdoodle_sketcher, name='chemdoodle_sketcher'),
    path('get_doi_info/<path:doi>/', utility_views.get_doi_info, name='get_doi_info'),

    path('get_graph_info/', graph_views.get_graph_info, name='get_graph_info'),

] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
