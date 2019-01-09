# Shadowfax website

The [Shadowfax website](http://www.dwarfs.ugent.be/shadowfax) makes use of
Jekyll, a Ruby based framework that generates a static HTML website from a
number of source files. The advantage is that (a) you don't need to know a lot
of HTML or anything to expand the website, as making a new post is as easy as
writing a Wikipedia page, (b) the website can run on about any server, since the
server just has to send the HTML pages to the client browser, and all fancy
stuff is done on the client side using Javascript and HTML5. This is different
for most other popular frameworks, which require a database living on the
server.

The particular framework used for this website is based on
<https://github.com/scotch-io/scotch-io.github.io> and has been adapted to our
liking. The most important differences are:

- a higher navigation bar showing the Shadowfax logo
- an image carousel on the front page
- an automated *Gallery* page with a fancy pop up window to show images and
video full size
- automatic outlining of test problem files on the *Test problems* page

The code should be quite self-explanatory. To set up the required Jekyll
development environment to make changes to the site, see
<https://scotch.io/tutorials/getting-started-with-jekyll-plus-a-free-bootstrap-3-starter-theme>.

To run the site locally during development, adapt the `url` line in
`global.yml`. Make sure to restore it to the public `url` before you publish the
site (and before running the final `jekyll build`). To publish the site, just
copy the entire contents of `_site/` to
`//files.ugent.be/USERNAME/www/shares/dwarfs/WWW/shadowfax/`.
