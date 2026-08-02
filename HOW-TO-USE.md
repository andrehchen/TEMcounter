# How to measure particle sizes with TEMcounter

You do not need to know anything about programming. You will copy two commands and
change one thing in each.

**Find your section below — Windows or Mac — and ignore the other one.**

---

# PART 1 — Setup

You only ever do this once. Budget 20 minutes, mostly waiting for downloads.

## Windows setup

**1. Install Python.** Go to <https://www.python.org/downloads/> and click the big
yellow download button. Run the installer.

> On the very first screen there is a checkbox at the bottom: **"Add python.exe to
> PATH"**. **Tick it before clicking Install.** If you miss it, nothing else will
> work and you will have to reinstall.

**2. Install Tesseract.** Go to
<https://github.com/UB-Mannheim/tesseract/wiki> and download the installer near the
top of the page. Run it and accept every default. This is what reads the "200 nm"
label off your images.

**3. Put the TEMcounter folder on your Desktop.** If you were sent a zip file,
right-click it and choose "Extract All".

**4. Open a terminal in that folder.** Open the TEMcounter folder in File Explorer.
Click once in the address bar at the top (where the folder path is), type
`powershell`, and press Enter. A blue window opens. This is the terminal, and it is
already pointed at the right folder.

**5. Paste these two commands.** Right-click pastes into PowerShell. Press Enter
after each, and wait for it to finish before the next.

```
python -m venv .venv
```

```
.venv\Scripts\pip install -r requirements.txt
```

The second one downloads a few hundred MB and prints a lot of text. That is normal.
When it stops and you get a blinking cursor back, setup is done.

## Mac setup

**1. Install Python.** Go to <https://www.python.org/downloads/> and click the big
yellow download button. Run the installer and accept the defaults.

**2. Open a terminal in the TEMcounter folder.** Press `Cmd + Space`, type
`Terminal`, press Enter. A window opens. Type `cd ` (with a space after it), then
**drag the TEMcounter folder from Finder into the Terminal window** — the path
fills itself in. Press Enter.

**3. Install Tesseract.** This reads the "200 nm" label off your images. Paste this
first and press Enter — it installs Homebrew, and will ask for your Mac password
(nothing appears as you type it, that is normal):

```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Then paste this:

```
brew install tesseract
```

**4. Paste these two commands**, pressing Enter after each and waiting for each to
finish:

```
python3 -m venv .venv
```

```
.venv/bin/pip install -r requirements.txt
```

The second downloads a few hundred MB and prints a lot of text. That is normal.

---

# PART 2 — Using it

Every time you want to measure images, first open the terminal in the TEMcounter
folder again, the same way you did in setup:

- **Windows:** open the TEMcounter folder, click the address bar, type `powershell`,
  Enter.
- **Mac:** open Terminal, type `cd ` and drag the TEMcounter folder in, Enter.

Then use one of the two commands below.

## Measuring ONE image

**Windows**

```
.venv\Scripts\python temcounter.py "PASTE FILE HERE" -n 150
```

**Mac**

```
.venv/bin/python temcounter.py "PASTE FILE HERE" -n 150
```

To fill in the file: delete `PASTE FILE HERE` (keep the quote marks), then

- **Windows:** hold **Shift**, right-click your `.tif` file, choose **"Copy as
  path"**, and paste it in. That already includes quote marks, so delete the ones
  in the command first.
- **Mac:** drag the `.tif` file from Finder straight into the Terminal window.

Press Enter.

## Measuring a WHOLE FOLDER

Exactly the same, but give it the folder instead of the file:

**Windows**

```
.venv\Scripts\python temcounter.py "PASTE FOLDER HERE" -n 150
```

**Mac**

```
.venv/bin/python temcounter.py "PASTE FOLDER HERE" -n 150
```

Fill it in the same way — Shift + right-click the folder and "Copy as path" on
Windows, or drag the folder into Terminal on a Mac.

`-n [integer]` is how many particles to measure per image. Change the number if you want
more or fewer.

## What happens next

It asks:

```
Name this run [260802_Summary1]:
```

**Just press Enter.** It names the run by today's date. Or type a name if you want
one, then press Enter.

Then it works for about 15 seconds per image and prints the answer:

```
diameter: 44.9 +/- 3.9 nm  (median 45.6, range 34.3-53.5, CV 8.6%, n=105)
```

That is your result. Everything is also saved in a `results` folder inside the
TEMcounter folder, in a subfolder named after the run:

| file | what it is |
|---|---|
| `<run name>.csv` | the summary — open it in Excel |
| `<run name>_particles.csv` | every single particle measured |
| `overlays` | pictures of your images with the measured particles circled |
| `histograms` | the size distribution graph |

**Always open one of the overlay pictures.** If the red circles sit neatly on the
particles, the numbers are good. That is the only check you need to do.

---

# If something goes wrong

**"python is not recognized"** (Windows) — Python was installed without ticking
"Add python.exe to PATH". Reinstall it and tick the box.

**"Tesseract OCR ... was not found"** — Step 2 of setup was skipped, or it went
somewhere unusual. Reinstall Tesseract and accept the default location.

**"No such file or directory"** — the file path is wrong. Do not type it by hand;
use "Copy as path" (Windows) or drag the file in (Mac).

**Red text mentioning the scale bar** — the tool could not read the scale bar in
that image. Check the image actually has one. If it does, send it to Andre.

**Nothing happens and the cursor just blinks** — it is working. Large images take
a while. Wait.

Anything else, send a screenshot of the whole terminal window to Andre.
