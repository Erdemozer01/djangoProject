import mimetypes
import os
from pathlib import Path
from django.http import HttpResponse
from django.shortcuts import render, reverse


def swiss_download(request):
    # Define Django project base directory
    BASE_DIR = Path(__file__).resolve().parent.parent
    # Define text file name
    filename = 'swiss_prot.txt'
    # Define the full file path
    filepath = os.path.join(BASE_DIR, 'files\\swiss_prot.txt')
    # Open the file for reading content
    path = open(filepath, 'r')
    # Set the mime type
    mime_type, _ = mimetypes.guess_type(filepath)
    # Set the return value of the HttpResponse
    response = HttpResponse(path, content_type=mime_type)
    # Set the HTTP header for sending to browser
    response['Content-Disposition'] = "attachment; filename=%s" % filename
    # Return the response value

    try:

        return response

    except FileNotFoundError:

        msg = "İndirmeye çalıştığınız dosya bulunamadı"

        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def fasta_download(request):
    # Define Django project base directory
    BASE_DIR = Path(__file__).resolve().parent.parent
    # Define text file name
    filename = 'file.fasta'
    # Define the full file path
    filepath = os.path.join(BASE_DIR, 'files\\file.fasta')
    # Open the file for reading content
    path = open(filepath, 'r')
    # Set the mime type
    mime_type, _ = mimetypes.guess_type(filepath)
    # Set the return value of the HttpResponse
    response = HttpResponse(path, content_type=mime_type)
    # Set the HTTP header for sending to browser
    response['Content-Disposition'] = "attachment; filename=%s" % filename
    # Return the response value

    try:

        return response

    except FileNotFoundError:

        msg = "İndirmeye çalıştığınız dosya bulunamadı"

        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def genbank_download(request):
    # Define Django project base directory
    BASE_DIR = Path(__file__).resolve().parent.parent
    # Define text file name
    filename = 'file.gbk'
    # Define the full file path
    filepath = os.path.join(BASE_DIR, 'files\\file.gbk')
    # Open the file for reading content
    path = open(filepath, 'r')
    # Set the mime type
    mime_type, _ = mimetypes.guess_type(filepath)
    # Set the return value of the HttpResponse
    response = HttpResponse(path, content_type=mime_type)
    # Set the HTTP header for sending to browser
    response['Content-Disposition'] = "attachment; filename=%s" % filename
    # Return the response value

    try:
        return response

    except FileNotFoundError:

        msg = "İndirmeye çalıştığınız dosya bulunamadı"

        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def hsp_download(request):
    try:
        # Define Django project base directory
        BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # Define text file name
        filename = 'hsp.txt'
        # Define the full file path
        filepath = BASE_DIR + '\\files\\' + filename
        # Open the file for reading content
        path = open(filepath, 'r')

        # Set the mime type
        mime_type, _ = mimetypes.guess_type(filepath)
        # Set the return value of the HttpResponse
        response = HttpResponse(path, content_type=mime_type)
        # Set the HTTP header for sending to browser
        response['Content-Disposition'] = "attachment; filename=%s" % filename
        # Return the response value
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        url = reverse \
            ("bioinformatic:xml_file")
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg, 'url': url})
    try:
        return response
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def entrez_download(request):
    try:
        # Define Django project base directory
        BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # Define text file name
        filename = 'entrez.txt'
        # Define the full file path
        filepath = BASE_DIR + '\\files\\' + filename
        # Open the file for reading content
        path = open(filepath, 'r')

        # Set the mime type
        mime_type, _ = mimetypes.guess_type(filepath)
        # Set the return value of the HttpResponse
        response = HttpResponse(path, content_type=mime_type)
        # Set the HTTP header for sending to browser
        response['Content-Disposition'] = "attachment; filename=%s" % filename
        # Return the response value
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        url = reverse("bioinformatic:entrez_file_search")
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg, 'url': url})
    try:
        return response
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def tree_download(request):
    try:
        # Define Django project base directory
        BASE_DIR = Path(__file__).resolve().parent.parent.parent
        # Define text file name
        filename = "{}_filogenetik_ağaç.jpg".format(request.user.username)
        # Define the full file path
        filepath = os.path.join(BASE_DIR, "media", "msa", "{}".format(request.user),
                                "{}_filogenetik_ağaç.jpg".format(request.user.username))
        # Open the file for reading content
        path = open(filepath, 'rb').read()
        # Set the mime type
        mime_type, _ = mimetypes.guess_type(filepath)
        # Set the return value of the HttpResponse
        response = HttpResponse(path, content_type=mime_type)
        # Set the HTTP header for sending to browser
        response['Content-Disposition'] = "attachment; filename=%s" % filename
        # Return the response value
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        url = reverse("bioinformatic:multiple_sequence_alignments")
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg, 'url': url})
    try:
        return response
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def clustal_alignment_download(request):
    try:
        # Define Django project base directory
        BASE_DIR = Path(__file__).resolve().parent.parent.parent
        # Define text file name
        from bioinformatic.models import MultipleSequenceAlignment
        get_file_type = MultipleSequenceAlignment.objects.all().filter(user=request.user.id).latest(
            'created').alignment_filetype
        # Define text file name
        filename = "{}_aligned.{}".format(request.user.username, get_file_type)
        # Define the full file path
        filepath = os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user.username),
                                "{}_aligned.{}".format(request.user.username, get_file_type))
        path = open(filepath, 'rb').read()
        # Set the mime type
        mime_type, _ = mimetypes.guess_type(filepath)
        # Set the return value of the HttpResponse
        response = HttpResponse(path, content_type=mime_type)
        # Set the HTTP header for sending to browser
        response['Content-Disposition'] = "attachment; filename=%s" % filename
        # Return the response value
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        url = reverse("bioinformatic:multiple_sequence_alignments")
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg, 'url': url})
    try:
        return response
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def maximum_likelihood_download(request):
    try:
        # Define Django project base directory
        BASE_DIR = Path(__file__).resolve().parent.parent.parent
        # Define text file name
        from bioinformatic.models import MultipleSequenceAlignment
        get_file_type = MultipleSequenceAlignment.objects.all().filter(user=request.user.id).latest(
            'created').palm_tools
        # Define text file name
        filename = "{}_result_{}.txt".format(request.user.username, get_file_type)
        # Define the full file path
        filepath = os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user.username),
                                "{}_result_{}.txt".format(request.user.username, get_file_type))
        path = open(filepath, 'rb').read()
        # Set the mime type
        mime_type, _ = mimetypes.guess_type(filepath)
        # Set the return value of the HttpResponse
        response = HttpResponse(path, content_type=mime_type)
        # Set the HTTP header for sending to browser
        response['Content-Disposition'] = "attachment; filename=%s" % filename
        # Return the response value
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        url = reverse("bioinformatic:maximum_likelihood")
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg, 'url': url})
    try:
        return response
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def PhyloXML_download(request):
    try:
        # Define Django project base directory
        BASE_DIR = Path(__file__).resolve().parent.parent.parent
        # Define text file name
        from bioinformatic.models import MultipleSequenceAlignment
        get_file_type = MultipleSequenceAlignment.objects.all().filter(user=request.user.id).latest(
            'created').palm_tools
        # Define text file name
        filename = "{}_tree.xml".format(request.user.username)
        # Define the full file path
        filepath = os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user.username),
                                "{}_tree.xml".format(request.user.username))
        path = open(filepath, 'rb').read()
        # Set the mime type
        mime_type, _ = mimetypes.guess_type(filepath)
        # Set the return value of the HttpResponse
        response = HttpResponse(path, content_type=mime_type)
        # Set the HTTP header for sending to browser
        response['Content-Disposition'] = "attachment; filename=%s" % filename
        # Return the response value
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        url = reverse("bioinformatic:maximum_likelihood")
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg, 'url': url})
    try:
        return response
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)

def clustalomega_alignment_download(request):
    try:
        # Define Django project base directory
        BASE_DIR = Path(__file__).resolve().parent.parent.parent
        # Define text file name
        from bioinformatic.models import MultipleSequenceAlignment
        get_file_type = MultipleSequenceAlignment.objects.all().filter(user=request.user.id).latest(
            'created').alignment_filetype
        # Define text file name
        filename = "{}_aligned.{}".format(request.user.username, get_file_type)
        # Define the full file path
        filepath = os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user.username),
                                "{}_aligned.{}".format(request.user.username, get_file_type))
        path = open(filepath, 'rb').read()
        # Set the mime type
        mime_type, _ = mimetypes.guess_type(filepath)
        # Set the return value of the HttpResponse
        response = HttpResponse(path, content_type=mime_type)
        # Set the HTTP header for sending to browser
        response['Content-Disposition'] = "attachment; filename=%s" % filename
        # Return the response value
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        url = reverse("bioinformatic:multiple_sequence_alignments")
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg, 'url': url})
    try:
        return response
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def clustal_stats_download(request):
    try:
        # Define Django project base directory
        BASE_DIR = Path(__file__).resolve().parent.parent.parent
        # Define text file name
        filename = "{}_stats.txt".format(request.user)
        # Define the full file path
        filepath = os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user),
                                "{}_stats.txt".format(request.user))
        # Open the file for reading content
        path = open(filepath, 'rb').read()
        # Set the mime type
        mime_type, _ = mimetypes.guess_type(filepath)
        # Set the return value of the HttpResponse
        response = HttpResponse(path, content_type=mime_type)
        # Set the HTTP header for sending to browser
        response['Content-Disposition'] = "attachment; filename=%s" % filename
        # Return the response value
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        url = reverse("bioinformatic:multiple_sequence_alignments")
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg, 'url': url})
    try:
        return response
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def clustal_scores_download(request):
    try:
        # Define Django project base directory
        BASE_DIR = Path(__file__).resolve().parent.parent.parent
        # Define text file name
        filename = "{}_scores.txt".format(request.user.username)
        # Define the full file path
        filepath = os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user),
                                "{}_scores.txt".format(request.user.username))
        # Open the file for reading content
        path = open(filepath, 'rb').read()
        # Set the mime type
        mime_type, _ = mimetypes.guess_type(filepath)
        # Set the return value of the HttpResponse
        response = HttpResponse(path, content_type=mime_type)
        # Set the HTTP header for sending to browser
        response['Content-Disposition'] = "attachment; filename=%s" % filename
        # Return the response value
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        url = reverse("bioinformatic:multiple_sequence_alignments")
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg, 'url': url})
    try:
        return response
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def global_alignments_download(request):
    # Define Django project base directory
    BASE_DIR = Path(__file__).resolve().parent.parent
    # Define text file name
    filename = 'global_alignment.txt'
    # Define the full file path
    filepath = os.path.join(BASE_DIR, 'files\\global_alignment.txt')
    # Open the file for reading content
    path = open(filepath, 'r')
    # Set the mime type
    mime_type, _ = mimetypes.guess_type(filepath)
    # Set the return value of the HttpResponse
    response = HttpResponse(path, content_type=mime_type)
    # Set the HTTP header for sending to browser
    response['Content-Disposition'] = "attachment; filename=%s" % filename
    # Return the response value
    try:
        return response
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def local_alignments_download(request):
    # Define Django project base directory
    BASE_DIR = Path(__file__).resolve().parent.parent

    # Define text file name
    filename = 'local_alignment.txt'
    # Define the full file path
    filepath = os.path.join(BASE_DIR, 'files\\local_alignment.txt')
    # Open the file for reading content
    path = open(filepath, 'r')
    # Set the mime type
    mime_type, _ = mimetypes.guess_type(filepath)
    # Set the return value of the HttpResponse
    response = HttpResponse(path, content_type=mime_type)
    # Set the HTTP header for sending to browser
    response['Content-Disposition'] = "attachment; filename=%s" % filename
    # Return the response value
    try:
        return response
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def muscle_aligned_download(request):
    # Define Django project base directory
    BASE_DIR = Path(__file__).resolve().parent.parent.parent
    from bioinformatic.models import MultipleSequenceAlignment
    get_file_type = MultipleSequenceAlignment.objects.all().filter(user=request.user.id).latest(
        'created').alignment_filetype
    # Define text file name
    filename = "{}_aligned.{}".format(request.user.username, get_file_type)
    # Define the full file path
    filepath = os.path.join(BASE_DIR, "media", 'msa', '{}'.format(request.user.username),
                            "{}_aligned.{}".format(request.user.username, get_file_type))
    # Open the file for reading content
    path = open(filepath, 'r')
    # Set the mime type
    mime_type, _ = mimetypes.guess_type(filepath)
    # Set the return value of the HttpResponse
    response = HttpResponse(path, content_type=mime_type)
    # Set the HTTP header for sending to browser
    response['Content-Disposition'] = "attachment; filename=%s" % filename
    # Return the response value
    try:
        return response
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        return render(request, 'bioinformatic/fasta/notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)
