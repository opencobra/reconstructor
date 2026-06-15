"""Views for collecting user feedback.

Feedback submitted through the UI is:
  1. Stored locally as a :class:`~reactions.models.Feedback` object so admins
     can browse it in the Django admin.
  2. Mirrored as an issue on the configured GitHub repository
     (``settings.GITHUB_FEEDBACK_REPO``) when a token is available.

GitHub mirroring is best-effort: if it fails (no token, network error, bad
credentials) the feedback is still saved and the API responds successfully,
recording the failure on the feedback object for later inspection.
"""

import json
import logging

import requests
from django.conf import settings
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt

from reactions.models import Feedback, User

logger = logging.getLogger(__name__)

# Valid feedback types, derived from the model so the two stay in sync.
VALID_FEEDBACK_TYPES = {choice[0] for choice in Feedback.FEEDBACK_TYPE_CHOICES}


def _create_github_issue(feedback):
    """Create a GitHub issue mirroring ``feedback``.

    Returns a tuple ``(issue_url, issue_number)`` on success. Raises an
    exception (with a human-readable message) on failure so the caller can
    record it.
    """
    token = getattr(settings, 'GITHUB_FEEDBACK_TOKEN', '')
    repo = getattr(settings, 'GITHUB_FEEDBACK_REPO', '')
    if not token:
        raise RuntimeError('GITHUB_FEEDBACK_TOKEN is not configured')
    if not repo:
        raise RuntimeError('GITHUB_FEEDBACK_REPO is not configured')

    type_label = feedback.get_feedback_type_display()
    submitter = feedback.user.name if feedback.user else 'Anonymous'

    title = f"[{type_label}] {feedback.text.strip().splitlines()[0][:80]}"
    body = (
        f"{feedback.text.strip()}\n\n"
        f"---\n"
        f"*Submitted via the Reconstructor feedback form.*\n"
        f"- **Type:** {type_label}\n"
        f"- **Submitted by:** {submitter}\n"
        f"- **Local feedback id:** {feedback.id}\n"
    )

    url = f"https://api.github.com/repos/{repo}/issues"
    headers = {
        'Authorization': f'Bearer {token}',
        'Accept': 'application/vnd.github+json',
        'X-GitHub-Api-Version': '2022-11-28',
    }
    payload = {'title': title, 'body': body, 'labels': [feedback.feedback_type]}

    response = requests.post(url, headers=headers, json=payload, timeout=15)
    if response.status_code not in (200, 201):
        raise RuntimeError(
            f'GitHub API returned {response.status_code}: {response.text[:300]}')

    data = response.json()
    return data.get('html_url'), data.get('number')


@csrf_exempt
def submit_feedback(request):
    """Persist user feedback and mirror it to GitHub.

    Expects a JSON body with:
        - ``feedback_type`` (str): one of the model's choice keys.
        - ``text`` (str): the feedback content.
        - ``user_id`` (int, optional): the submitting user's id.

    Returns:
        JsonResponse describing whether the feedback was saved and whether the
        GitHub issue was created.
    """
    if request.method != 'POST':
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid request method'}, status=405)

    try:
        data = json.loads(request.body or '{}')
    except json.JSONDecodeError:
        return JsonResponse(
            {'status': 'error', 'message': 'Invalid JSON body'}, status=400)

    feedback_type = (data.get('feedback_type') or 'other').strip()
    text = (data.get('text') or '').strip()
    user_id = data.get('user_id')

    if feedback_type not in VALID_FEEDBACK_TYPES:
        feedback_type = 'other'
    if not text:
        return JsonResponse(
            {'status': 'error', 'message': 'Feedback text is required'},
            status=400)

    user = None
    if user_id:
        user = User.objects.filter(pk=user_id).first()

    feedback = Feedback.objects.create(
        feedback_type=feedback_type, text=text, user=user)

    # Best-effort GitHub mirroring.
    github_created = False
    try:
        issue_url, issue_number = _create_github_issue(feedback)
        feedback.github_issue_url = issue_url
        feedback.github_issue_number = issue_number
        feedback.save(update_fields=['github_issue_url', 'github_issue_number'])
        github_created = True
    except Exception as exc:  # noqa: BLE001 - mirror failures must not break save
        logger.warning('Failed to create GitHub issue for feedback %s: %s',
                       feedback.id, exc)
        feedback.github_sync_error = str(exc)
        feedback.save(update_fields=['github_sync_error'])

    return JsonResponse({
        'status': 'success',
        'message': 'Thank you for your feedback!',
        'feedback_id': feedback.id,
        'github_issue_created': github_created,
        'github_issue_url': feedback.github_issue_url,
    })
