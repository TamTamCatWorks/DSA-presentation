## Get Quarto

Detailed instructions at [quarto docs](https://quarto.org/docs/get-started/).

## Build

### Main Slides

```
quarto render main.qmd
```

Or preview live while editing:

```
quarto preview main.qmd
```

### Individual Sections

```
cd slides
quarto render {section}.qmd
```

with `{section}` your desired section

## View Slides

Open the built `.html` file(s) in the browser and it should display the slides. 

```
<!-- Linux -->
xdg-open _site/index.html 
<!-- macOS  -->
open _site/index.html 
<!-- Windows -->
start _site/index.html 
```

Or just drag the `.html` into your browser.

## GitHub Pages

Every push to `main` builds and publishes the combined presentation and each
standalone section with GitHub Actions. The published site is available at:

https://tamtamcatworks.github.io/DSA-presentation/

The combined presentation is the main page. Standalone sections are published
under `/sections/`, for example:

https://tamtamcatworks.github.io/DSA-presentation/sections/02-hand-trace/

Each workflow run also provides a downloadable `slides-site` artifact from the
Actions run summary. Generated files are kept out of the repository.

## Remarks

- I gitingore the built files since I don't want to clutter the repo, so you have to build them yourself to see the work.
- Usually I use the built slides as a starting point to load in your desired editor (like powerpoint or google slides). I add images there, then convert it back to markdown with [pptx2md](https://github.com/ssine/pptx2md). Repeat the process until it's good.

### Pasting Images

- One thing that VSCode can't do is paste images directly into `.qmd` files, so it can be annoying when editing.
- To fix this you can either:
    - Enter Visual Editing Mode
    - Temporarily change the file extension into `.md`

Personally I just use RStudio. 