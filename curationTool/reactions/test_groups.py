"""
Backend tests for the multi-reaction tabbed workspace (001-multi-reaction-tabbed).

Covers the persistence-critical group endpoints:
- group-contents returns ordered members (SC-001 / T015)
- create / rename / delete group; delete keeps member reactions (FR-009 / T021)
- add persists, remove keeps the saved reaction, duplicate add guarded
  (SC-003 / FR-011 / FR-012 / T025)
"""
import json

from django.test import TestCase
from django.urls import reverse

from reactions.models import (
    Flag,
    Reaction,
    ReactionGroup,
    ReactionGroupMembership,
    User,
)


def make_reaction(name):
    return Reaction.objects.create(
        substrates='a', products='b', short_name=name)


class GroupBackendTests(TestCase):
    def setUp(self):
        self.user = User.objects.create(name='curator')
        self.r1 = make_reaction('R1')
        self.r2 = make_reaction('R2')
        self.r3 = make_reaction('R3')

    # --- helpers -------------------------------------------------------
    def post(self, name, payload):
        return self.client.post(
            reverse(name), data=json.dumps(payload),
            content_type='application/json')

    def add(self, group, reactions):
        return self.post('add_to_group', {
            'userID': self.user.pk, 'groupId': group.pk,
            'reactionIds': [r.pk for r in reactions]})

    # --- T015 / SC-001 -------------------------------------------------
    def test_contents_returns_ordered_members(self):
        group = self.user.reaction_groups.create(name='G', is_active=True)
        self.add(group, [self.r1, self.r2, self.r3])
        resp = self.client.get(
            reverse('group_contents'),
            {'userID': self.user.pk, 'groupId': group.pk})
        body = resp.json()
        ids = [c['id'] for c in body['reactions']]
        self.assertEqual(ids, [self.r1.pk, self.r2.pk, self.r3.pk])
        self.assertEqual(body['group']['count'], 3)

    def test_default_group_seeded_when_none(self):
        # FR-008 / SC-005: a user with no groups still gets one.
        resp = self.client.get(
            reverse('group_contents'), {'userID': self.user.pk})
        self.assertEqual(resp.json()['status'], 'success')
        self.assertEqual(self.user.reaction_groups.count(), 1)
        self.assertTrue(self.user.reaction_groups.first().is_active)

    # --- T021 / FR-009 -------------------------------------------------
    def test_create_rename_delete_group(self):
        created = self.post(
            'create_group', {'userID': self.user.pk, 'name': 'Working set'})
        gid = created.json()['group']['id']
        self.assertEqual(ReactionGroup.objects.get(pk=gid).name, 'Working set')

        self.post('rename_group',
                  {'userID': self.user.pk, 'groupId': gid, 'name': 'Renamed'})
        self.assertEqual(ReactionGroup.objects.get(pk=gid).name, 'Renamed')

        self.post('delete_group', {'userID': self.user.pk, 'groupId': gid})
        self.assertFalse(ReactionGroup.objects.filter(pk=gid).exists())

    def test_delete_group_keeps_member_reactions(self):
        group = self.user.reaction_groups.create(name='G', is_active=True)
        self.add(group, [self.r1, self.r2])
        self.post('delete_group', {'userID': self.user.pk, 'groupId': group.pk})
        # FR-009: reactions survive group deletion.
        self.assertTrue(Reaction.objects.filter(pk=self.r1.pk).exists())
        self.assertTrue(Reaction.objects.filter(pk=self.r2.pk).exists())

    def test_always_at_least_one_group_after_delete(self):
        group = self.user.reaction_groups.create(name='only', is_active=True)
        self.post('delete_group', {'userID': self.user.pk, 'groupId': group.pk})
        self.assertGreaterEqual(self.user.reaction_groups.count(), 1)  # FR-008

    # --- T025 / SC-003 / FR-011 / FR-012 -------------------------------
    def test_add_persists_membership(self):
        group = self.user.reaction_groups.create(name='G', is_active=True)
        self.add(group, [self.r1])
        self.assertTrue(ReactionGroupMembership.objects.filter(
            group=group, reaction=self.r1).exists())

    def test_remove_keeps_saved_reaction(self):
        group = self.user.reaction_groups.create(name='G', is_active=True)
        self.add(group, [self.r1])
        self.post('remove_from_group',
                  {'userID': self.user.pk, 'groupId': group.pk,
                   'reactionId': self.r1.pk})
        self.assertFalse(ReactionGroupMembership.objects.filter(
            group=group, reaction=self.r1).exists())
        # FR-011: removal does not delete the underlying reaction.
        self.assertTrue(Reaction.objects.filter(pk=self.r1.pk).exists())

    def test_duplicate_add_prevented(self):
        group = self.user.reaction_groups.create(name='G', is_active=True)
        self.add(group, [self.r1])
        resp = self.add(group, [self.r1])
        self.assertIn(self.r1.pk, resp.json()['already_present'])
        # FR-012: still exactly one membership row.
        self.assertEqual(ReactionGroupMembership.objects.filter(
            group=group, reaction=self.r1).count(), 1)

    def test_saved_reactions_for_group_includes_filter_metadata(self):
        group = self.user.reaction_groups.create(name='G', is_active=True)
        flag = Flag.objects.create(
            user=self.user, name_flag='Review', color='#ff0000')
        rxn = Reaction.objects.create(
            substrates='atp[c]', products='adp[c]', short_name='ATPase',
            rxn_formula='atp[c] -> adp[c]', subsystem='Energy',
            direction='forward')
        rxn.flags.add(flag)
        self.user.saved_reactions.add(rxn)
        self.add(group, [rxn])

        resp = self.client.get(
            reverse('saved_reactions_for_group'),
            {'userID': self.user.pk, 'groupId': group.pk})
        card = resp.json()['reactions'][0]
        self.assertEqual(card['subsystem'], 'Energy')
        self.assertEqual(card['substrates'], 'atp[c]')
        self.assertEqual(card['products'], 'adp[c]')
        self.assertEqual(card['direction'], 'forward')
        self.assertTrue(card['in_group'])
        self.assertEqual(card['flags'][0]['name'], 'Review')
        self.assertEqual(card['flags'][0]['color'], '#ff0000')

    def test_anonymous_contents_is_safe(self):
        # FR-018: no logged-in user → no error, empty payload.
        resp = self.client.get(reverse('group_contents'))
        self.assertEqual(resp.json()['status'], 'anonymous')

    # --- T029 / SC-004 / FR-014 / FR-015 -------------------------------
    def test_clone_yields_distinct_reaction_original_unchanged(self):
        """
        Cloning copies a reaction's fields into a new row; editing + saving the
        clone (different metabolite + GPR) must not alter the original, and the
        two must be distinct persisted reactions.
        """
        original = Reaction.objects.create(
            substrates='atp[c]', products='adp[c]', short_name='ORIG',
            rxn_formula='atp[c] -> adp[c]', gene_info=[{'gene': 'G1'}],
            subsystem='Transport')
        # Snapshot the original's persisted fields (SC-004: byte-for-byte).
        before = {
            'substrates': original.substrates,
            'products': original.products,
            'rxn_formula': original.rxn_formula,
            'gene_info': original.gene_info,
        }

        # Clone the way clone_reaction_view does: copy row, new pk.
        clone = Reaction.objects.get(pk=original.pk)
        clone.pk = None
        clone.short_name = 'ORIG (clone)'
        clone.save()
        # Edit the clone: swap a metabolite + change the GPR.
        clone.products = 'amp[c]'
        clone.gene_info = [{'gene': 'G2'}]
        clone.save()
        self.user.saved_reactions.add(clone)

        group = self.user.reaction_groups.create(name='G', is_active=True)
        self.add(group, [clone])

        # Distinct rows, both persisted.
        self.assertNotEqual(clone.pk, original.pk)
        self.assertTrue(Reaction.objects.filter(pk=clone.pk).exists())

        # Original byte-for-byte unchanged (FR-014).
        original.refresh_from_db()
        self.assertEqual(original.substrates, before['substrates'])
        self.assertEqual(original.products, before['products'])
        self.assertEqual(original.rxn_formula, before['rxn_formula'])
        self.assertEqual(original.gene_info, before['gene_info'])

        # The clone carries the edits.
        self.assertEqual(Reaction.objects.get(pk=clone.pk).products, 'amp[c]')

    # --- Card formula cleaning (fix #2) --------------------------------
    def test_contents_cleans_list_wrapped_formula(self):
        group = self.user.reaction_groups.create(name='G', is_active=True)
        rxn = Reaction.objects.create(
            substrates='a', products='b', short_name='F',
            rxn_formula='["atp[c] -> adp[c] + pi[c]"]')
        self.add(group, [rxn])
        resp = self.client.get(
            reverse('group_contents'),
            {'userID': self.user.pk, 'groupId': group.pk})
        card = resp.json()['reactions'][0]
        # The stored list-wrapper is stripped to a readable formula.
        self.assertEqual(card['rxn_formula'], 'atp[c] -> adp[c] + pi[c]')

    # --- Discard clone (fix #3 / FR-015) -------------------------------
    def test_discard_clone_deletes_row_and_membership(self):
        clone = make_reaction('CLONE')
        self.user.saved_reactions.add(clone)
        group = self.user.reaction_groups.create(name='G', is_active=True)
        self.add(group, [clone])
        self.post('discard_clone',
                  {'userID': self.user.pk, 'reactionId': clone.pk})
        self.assertFalse(Reaction.objects.filter(pk=clone.pk).exists())
        self.assertFalse(ReactionGroupMembership.objects.filter(
            reaction_id=clone.pk).exists())

    def test_discard_refuses_reaction_not_owned(self):
        # Safety: cannot delete a reaction not in the user's saved set.
        other = make_reaction('NOT_MINE')
        resp = self.post('discard_clone',
                         {'userID': self.user.pk, 'reactionId': other.pk})
        self.assertEqual(resp.status_code, 403)
        self.assertTrue(Reaction.objects.filter(pk=other.pk).exists())

    # --- Targeted gene-info session clear (fix #3) ---------------------
    def test_clear_gene_info_session_keeps_login(self):
        session = self.client.session
        session['userID'] = self.user.pk
        session['gene_info'] = [{'info': 'GENE1'}]
        session.save()
        resp = self.client.post(reverse('clear_gene_info_session'))
        self.assertEqual(resp.json()['status'], 'success')
        self.assertNotIn('gene_info', self.client.session)
        self.assertEqual(self.client.session['userID'], self.user.pk)
