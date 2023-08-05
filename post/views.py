from django.views.generic import ListView, DetailView, CreateView, UpdateView, DeleteView
from .models import Posts
from django.urls import reverse_lazy
from .models import Category
from django.shortcuts import get_object_or_404
from .forms import PostForm



class CategoriesView(ListView):
    template_name = 'blog/pages/categories.html'
    model = Category
    paginate_by = 10


class CategoryView(ListView):
    template_name = 'blog/pages/category.html'
    model = Posts
    paginate_by = 3

    def get_queryset(self):
        category = get_object_or_404(Category, slug=self.kwargs['slug'])
        return category.post.all().order_by('-id')


class PostDetailView(DetailView):
    template_name = 'blog/pages/detail.html'
    model = Posts
    queryset = Posts.objects.all()
    context_object_name = "post"


class AddPostView(CreateView):
    template_name = 'pages/add_post.html'
    model = Posts
    form_class = PostForm
    success_url = reverse_lazy('blog:home')


class AddCategoryView(CreateView):
    template_name = 'pages/add_category.html'
    model = Category
    fields = '__all__'
    success_url = reverse_lazy('post:add_post')


class UpdatePostView(UpdateView):
    template_name = 'pages/update_post.html'
    model = Posts
    form_class = PostForm
    success_url = reverse_lazy('blog:home')


class DeletePostView(DeleteView):
    template_name = 'pages/delete_post.html'
    model = Posts
    success_url = reverse_lazy('blog:home')
