# Publish Theta Atlas at kdwong.github.io

GitHub Pages serves static HTML, CSS, and JavaScript, but this calculator also
needs Python. The deployment therefore has two parts:

1. `kdwong.github.io` serves the files in `web/`.
2. A Python web service serves `POST /api/calculate`.

The browser performs the same certified calculation as the local version; it
does not use a precomputed answer table.

## 1. Put this calculator in a repository you control

This working tree comes from Jia-Jun Ma's
[`unipotentrepn`](https://github.com/jiajunma/unipotentrepn) repository. Create
a fork or another repository under the `kdwong` account, commit the calculator
files, and push them there. Do not try to deploy uncommitted local files.

The backend needs at least:

- `theta2_pbp.py`
- `theta2_pbp_web.py`
- `standalone.py`
- `combunipotent/`
- `web/`
- `requirements.txt`
- `requirements-local.txt`
- `render.yaml`

## 2. Deploy the Python API

One straightforward option is Render:

1. Sign in to Render with GitHub.
2. Create a new Blueprint and select the repository from step 1.
3. Render will read `render.yaml`, install the Python requirements, start the
   server, and check `/health`.
4. Copy the resulting address, such as
   `https://theta-atlas-api.onrender.com`.
5. Verify that opening
   `https://theta-atlas-api.onrender.com/health` shows `{"status":"ok"}`.

The server already accepts browser requests from `https://kdwong.github.io`.
To authorize another site, set the service environment variable
`THETA_PBP_ALLOWED_ORIGINS` to a comma-separated list of exact origins.

## 3. Add the interface to kdwong.github.io

Clone `https://github.com/kdwong/kdwong.github.io`, create a directory named
`theta-atlas`, and copy these three files into it:

- `web/index.html`
- `web/app.js`
- `web/styles.css`

In the copied `theta-atlas/index.html`, replace the empty API setting

```html
<meta name="theta-api-url" content="" />
```

with the deployed endpoint:

```html
<meta
  name="theta-api-url"
  content="https://theta-atlas-api.onrender.com/api/calculate"
/>
```

Commit and push the `theta-atlas` directory to the `main` branch. Because the
existing personal site is published from that repository, the calculator will
then be available at:

`https://kdwong.github.io/theta-atlas/`

Optionally add this link to the main `index.html` of the personal site:

```html
<a href="theta-atlas/">Theta Atlas: painted-bipartition calculator</a>
```

## 4. Final check

Open the public calculator in a private browser window, calculate `(4)`, and
confirm that:

- the final form is `O(2,3)`;
- both targets `(00|000)` and `(11|000)` appear, in that order;
- the three source links at the bottom open correctly.

The first request can take longer if a free backend has gone idle. A paid
always-on instance removes that startup delay.
