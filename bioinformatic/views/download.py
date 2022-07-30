import mimetypes
import os
from pathlib import Path
from django.http import HttpResponse
from django.shortcuts import render


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
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # Define text file name
    filename = 'file.gbk'
    # Define the full file path
    filepath = "media\\genbank.txt"
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
        return render(request, 'laboratory/biyoinformatik/dosya işlemleri/Writing Sequence Files/download.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def stockholm_download(request):
    # Define Django project base directory
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # Define text file name
    filename = 'file.sth'
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
    try:
        return response
    except FileNotFoundError:
        msg = "İndirmeye çalıştığınız dosya bulunamadı"
        return render(request, 'laboratory/biyoinformatik/dosya işlemleri/Writing Sequence Files/download.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def global_alignments_download(request):
    # Define Django project base directory
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # Define text file name
    filename = 'global.txt'
    # Define the full file path
    filepath = "media\\global.txt"
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
        return render(request, 'notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)


def local_alignments_download(request):
    # Define Django project base directory
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # Define text file name
    filename = 'local.txt'
    # Define the full file path
    filepath = "media\\local.txt"
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
        return render(request, 'notfound.html',
                      {"msg": msg})
    finally:
        os.remove(filepath)
