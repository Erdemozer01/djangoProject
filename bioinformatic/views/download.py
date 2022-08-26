import mimetypes
import os
from pathlib import Path
from django.http import HttpResponse
from django.shortcuts import render, reverse
from django.conf import settings


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
        url = reverse\
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
        # Define text file name
        filename = "tree.jpg"
        # Define the full file path
        filepath = os.path.join(settings.MEDIA_ROOT, "tree.jpg")
        # Open the file for reading content
        path = open(filepath, 'rb')
        # Set the mime type
        mime_type, _ = mimetypes.guess_type(filepath)
        # Set the return value of the HttpResponse
        response = HttpResponse(path, content_type=mime_type)
        # Set the HTTP header for sending to browser
        response['Content-Disposition'] = "attachment; filename=%s" % filename
        # Return the response value
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        url = reverse("bioinformatic:filogenetik_agac")
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
