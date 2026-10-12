## Get Quarto

Detailed instructions at ![quarto docs](https://quarto.org/docs/get-started/).

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

## Remarks

- I gitingore the built files since I don't want to clutter the repo, so you have to build them yourself to see the work.
- Usually I use the built slides as a starting point to load in your desired editor (like powerpoint or google slides). I add images there, then convert it back to markdown with ![pptx2md](https://github.com/ssine/pptx2md). Repeat the process until it's good.
