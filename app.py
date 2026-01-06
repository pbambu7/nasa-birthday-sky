# app.py
# NASA Birthday Sky — Streamlit app
# Requirements (requirements.txt):
# streamlit
# requests
# pillow
# astropy
# geopy

import io
import math
import random
import datetime as dt
from datetime import date, time as dtime

import requests
import streamlit as st
from PIL import Image, ImageDraw, ImageFont

from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, get_body
from astropy import units as u
from geopy.geocoders import Nominatim


# ---------------- PAGE CONFIG ----------------
st.set_page_config(
    page_title="NASA Birthday Sky",
    layout="centered",
)

st.title("🌌 NASA Birthday Sky")
st.caption("What was above you when you were born? Built with NASA open data + real sky calculations.")


# ---------------- SECRETS ----------------
NASA_API_KEY = st.secrets.get("NASA_API_KEY") or st.secrets.get("NASA_KEY")  # allow either name
if not NASA_API_KEY:
    st.error("NASA API key not found. Add it in Streamlit → App settings → Secrets:\n\nNASA_API_KEY = \"YOUR_KEY\"")
    st.stop()

# Optional: set this in Secrets so share links are perfect.
# APP_BASE_URL = "https://your-app-name.streamlit.app"
APP_BASE_URL = st.secrets.get("APP_BASE_URL", "").strip()


# ---------------- HELPERS ----------------
def _safe_str(x) -> str:
    return "" if x is None else str(x).strip()


def _parse_date(s: str):
    try:
        return dt.date.fromisoformat(s)
    except Exception:
        return None


def _parse_float(s: str):
    try:
        return float(s)
    except Exception:
        return None


def _clamp(n, lo, hi):
    return max(lo, min(hi, n))


# ---------------- URL PARAMS (share links) ----------------
# Using experimental_* for broad Streamlit compatibility
qp = st.experimental_get_query_params()
qp_date = _parse_date(_safe_str(qp.get("date", [None])[0]))
qp_city = _safe_str(qp.get("city", [""])[0])
qp_mode = _safe_str(qp.get("mode", ["midnight"])[0]).lower()  # midnight | noon | exact
qp_time = _safe_str(qp.get("time", ["00:00"])[0])  # HH:MM
qp_autorun = _safe_str(qp.get("autorun", ["0"])[0])  # 1 to autorun


# ---------------- UI INPUTS ----------------
col1, col2 = st.columns(2)

with col1:
    birth_date = st.date_input(
        "Your birth date",
        value=qp_date or date(2000, 1, 1),
        max_value=date.today(),
    )

with col2:
    birth_city = st.text_input(
        "Birth city",
        value=qp_city,
        placeholder="Salt Lake City",
    )

moment_mode = st.selectbox(
    "Moment to calculate (local time)",
    options=["Midnight (00:00)", "Noon (12:00)", "Exact time (enter below)"],
    index=0 if qp_mode == "midnight" else (1 if qp_mode == "noon" else 2),
)

exact_time_default = qp_time if qp_time else "00:00"
exact_time_str = st.text_input(
    "Exact birth time (HH:MM) — optional",
    value=exact_time_default,
    disabled=(moment_mode != "Exact time (enter below)"),
    help="If you don’t know the time, use Midnight or Noon.",
)

st.divider()


# ---------------- NETWORK: NASA APOD (robust) ----------------
@st.cache_data(show_spinner=False, ttl=60 * 60)  # cache for 1 hour
def fetch_apod_apodata(birth_date_iso: str, api_key: str):
    """
    Returns a dict (APOD JSON) or None if unavailable.
    Handles timeouts + retries + outage gracefully.
    """
    url = "https://api.nasa.gov/planetary/apod"
    params = {"date": birth_date_iso, "api_key": api_key}

    # NASA APOD is sometimes degraded; do small retries with backoff.
    timeouts = (4, 12)  # (connect, read)
    for attempt in range(4):
        try:
            r = requests.get(url, params=params, timeout=timeouts)
            r.raise_for_status()
            data = r.json()
            # APOD sometimes returns video; we handle that later.
            return data
        except requests.exceptions.Timeout:
            # exponential-ish backoff
            delay = 0.6 * (2 ** attempt)
            try:
                import time as _time
                _time.sleep(delay)
            except Exception:
                pass
        except Exception:
            # any other request error or JSON error
            return None
    return None


# ---------------- GEO: City -> coordinates (robust) ----------------
@st.cache_data(show_spinner=False, ttl=60 * 60 * 24)  # cache 24 hours
def geocode_city(city: str):
    """
    Returns (lat, lon, display_name) or None.
    """
    city = city.strip()
    if not city:
        return None

    try:
        geolocator = Nominatim(
            user_agent="nasa-birthday-sky",
            timeout=6,
        )
        loc = geolocator.geocode(city)
        if not loc:
            return None
        return float(loc.latitude), float(loc.longitude), _safe_str(loc.address)
    except Exception:
        return None


def resolve_birth_datetime(birth_date: date, mode_label: str, exact_time_str: str):
    """
    Returns a naive datetime (date + chosen time). We compute astro positions using astropy Time.
    NOTE: Without timezone, this is an approximation (still fun + consistent).
    """
    if mode_label.startswith("Midnight"):
        t = dtime(0, 0)
        mode = "midnight"
    elif mode_label.startswith("Noon"):
        t = dtime(12, 0)
        mode = "noon"
    else:
        mode = "exact"
        hh, mm = 0, 0
        try:
            parts = exact_time_str.strip().split(":")
            if len(parts) >= 2:
                hh = int(parts[0])
                mm = int(parts[1])
            hh = _clamp(hh, 0, 23)
            mm = _clamp(mm, 0, 59)
        except Exception:
            hh, mm = 0, 0
        t = dtime(hh, mm)

    return dt.datetime.combine(birth_date, t), mode


# ---------------- SKY CALC: visible planets ----------------
@st.cache_data(show_spinner=False, ttl=60 * 60)
def get_visible_planets(lat: float, lon: float, when_iso: str):
    """
    Returns list of dicts: [{"name":"Jupiter","alt_deg":..., "az_deg":...}, ...]
    Visible = altitude > 0 degrees.
    """
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)
    t = Time(when_iso)  # ISO string
    frame = AltAz(obstime=t, location=location)

    planets = ["venus", "mars", "jupiter"]  # keep it simple + fast
    out = []
    for p in planets:
        try:
            body = get_body(p, t).transform_to(frame)
            alt = body.alt.to(u.deg).value
            az = body.az.to(u.deg).value
            if alt > 0:
                out.append(
                    {"name": p.capitalize(), "alt_deg": float(alt), "az_deg": float(az)}
                )
        except Exception:
            continue
    # sort: highest first
    out.sort(key=lambda x: x["alt_deg"], reverse=True)
    return out


# ---------------- VISUALS: fallback starfield ----------------
def generate_starfield_ig(size=1080, seed_text=""):
    """
    Creates a shareable starfield PNG (PIL Image) deterministic from seed_text.
    """
    rnd = random.Random(seed_text)
    img = Image.new("RGB", (size, size), (0, 0, 0))
    draw = ImageDraw.Draw(img)

    # stars
    n = 1800
    for _ in range(n):
        x = rnd.randint(0, size - 1)
        y = rnd.randint(0, size - 1)
        # brightness distribution: many dim, few bright
        b = int(30 + (rnd.random() ** 0.35) * 225)
        r = 1 if rnd.random() < 0.92 else 2
        draw.ellipse((x - r, y - r, x + r, y + r), fill=(b, b, b))

    # subtle nebula-ish haze (very light)
    for _ in range(18):
        x = rnd.randint(0, size)
        y = rnd.randint(0, size)
        rad = rnd.randint(120, 280)
        alpha = rnd.randint(8, 16)
        haze = Image.new("RGBA", (size, size), (0, 0, 0, 0))
        hz = ImageDraw.Draw(haze)
        hz.ellipse((x - rad, y - rad, x + rad, y + rad), fill=(60, 80, 140, alpha))
        img = Image.alpha_composite(img.convert("RGBA"), haze).convert("RGB")

    return img


def create_share_card(
    base_img: Image.Image,
    birth_date: date,
    city: str,
    mode: str,
    visible_planets,
    apod_title: str = "",
    footer_note: str = "",
):
    """
    Returns BytesIO PNG ready to download/share.
    """
    img = base_img.copy().convert("RGB").resize((1080, 1080))
    draw = ImageDraw.Draw(img)

    overlay_h = 320
    draw.rectangle([(0, 1080 - overlay_h), (1080, 1080)], fill=(0, 0, 0))

    # Try a simple default font; Streamlit Cloud usually has a basic one.
    try:
        font_big = ImageFont.truetype("DejaVuSans.ttf", 44)
        font_mid = ImageFont.truetype("DejaVuSans.ttf", 32)
        font_small = ImageFont.truetype("DejaVuSans.ttf", 26)
    except Exception:
        font_big = font_mid = font_small = None

    planets_text = "None"
    if visible_planets:
        planets_text = ", ".join([p["name"] for p in visible_planets])

    header = "NASA Birthday Sky"
    line1 = f"Date: {birth_date.isoformat()}  •  City: {city}"
    line2 = f"Moment: {mode}"
    line3 = f"Planets above horizon: {planets_text}"

    y0 = 1080 - overlay_h + 34
    draw.text((48, y0), header, fill="white", font=font_big)
    draw.text((48, y0 + 70), line1, fill="white", font=font_mid)
    draw.text((48, y0 + 120), line2, fill=(200, 200, 200), font=font_small)
    draw.text((48, y0 + 160), line3, fill=(200, 200, 200), font=font_small)

    if apod_title:
        draw.text((48, y0 + 210), f"APOD: {apod_title}", fill=(200, 200, 200), font=font_small)

    if footer_note:
        draw.text((48, 1040), footer_note, fill=(160, 160, 160), font=font_small)

    buf = io.BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    return buf


def build_share_link(base_url: str, birth_date: date, city: str, mode: str, exact_time: str):
    """
    Returns a sharable URL with query params.
    If base_url not provided, returns just the query string.
    """
    q_city = city.strip().replace(" ", "%20")
    q = f"?date={birth_date.isoformat()}&city={q_city}&mode={mode}&time={exact_time}&autorun=1"
    if base_url:
        return base_url.rstrip("/") + "/" + q
    return q


# ---------------- ACTION ----------------
run_clicked = st.button("🚀 Run my birthday", use_container_width=True)
autorun = (qp_autorun == "1")
should_run = run_clicked or autorun

if should_run:
    if not birth_city.strip():
        st.warning("Please enter a birth city (example: Kinshasa, Salt Lake City).")
        st.stop()

    # Resolve city to coordinates
    with st.spinner("Finding your city coordinates..."):
        geo = geocode_city(birth_city)

    if not geo:
        st.error("Couldn’t find that city. Try adding country/state (e.g., “Kinshasa, Congo”).")
        st.stop()

    lat, lon, display_name = geo

    # Resolve moment
    when_dt, mode_key = resolve_birth_datetime(birth_date, moment_mode, exact_time_str)
    when_iso = when_dt.isoformat()

    # Get visible planets
    with st.spinner("Calculating planets above the horizon..."):
        visible_planets = get_visible_planets(lat, lon, when_iso)

    # Try NASA APOD (may be degraded)
    apod_data = None
    with st.spinner("Fetching NASA APOD (this can be temporarily degraded)..."):
        apod_data = fetch_apod_apodata(birth_date.isoformat(), NASA_API_KEY)

    apod_ok = False
    apod_img = None
    apod_title = ""
    apod_expl = ""
    apod_url = ""

    # Determine if APOD image is usable
    if apod_data and isinstance(apod_data, dict):
        apod_title = _safe_str(apod_data.get("title"))
        apod_expl = _safe_str(apod_data.get("explanation"))
        apod_url = _safe_str(apod_data.get("url"))
        media_type = _safe_str(apod_data.get("media_type", "image")).lower()

        if media_type == "image" and apod_url:
            try:
                # Download with retries; APOD image servers can also be slow
                img_bytes = None
                for attempt in range(3):
                    try:
                        rr = requests.get(apod_url, timeout=(4, 12))
                        rr.raise_for_status()
                        img_bytes = rr.content
                        break
                    except requests.exceptions.Timeout:
                        try:
                            import time as _time
                            _time.sleep(0.6 * (2 ** attempt))
                        except Exception:
                            pass
                    except Exception:
                        img_bytes = None
                        break

                if img_bytes:
                    apod_img = Image.open(io.BytesIO(img_bytes)).convert("RGB")
                    apod_ok = True
            except Exception:
                apod_ok = False

    # Render the visual section
    if apod_ok and apod_img is not None:
        st.success("✅ NASA APOD loaded successfully.")
        st.subheader(apod_title or "NASA Astronomy Picture of the Day")
        st.image(apod_img, use_column_width=True)
        if apod_expl:
            with st.expander("APOD explanation"):
                st.write(apod_expl)
    else:
        # Fallback starfield so the app never breaks during outages
        st.warning("🚨 NASA APOD is temporarily unavailable (NASA shows a service outage). Planetary calculations are still accurate.")
        seed_text = f"{birth_date.isoformat()}|{birth_city.strip().lower()}|{mode_key}|{exact_time_str}"
        fallback_img = generate_starfield_ig(seed_text=seed_text)
        st.image(fallback_img, use_column_width=True, caption="Generated starfield (fallback)")

    # Planets list
    st.subheader("🪐 Planets visible above the horizon")
    if visible_planets:
        for p in visible_planets:
            st.write(f"- **{p['name']}** — altitude **{p['alt_deg']:.1f}°**, azimuth **{p['az_deg']:.1f}°**")
    else:
        st.write("No major planets were above the horizon at that moment (for this simplified planet list).")

    st.caption(f"📍 Location used: {display_name} (lat {lat:.3f}, lon {lon:.3f})")

    # Share link
    st.subheader("🔗 Share this result")
    share_link = build_share_link(APP_BASE_URL, birth_date, birth_city, mode_key, exact_time_str)
    if APP_BASE_URL:
        st.code(share_link, language="text")
    else:
        st.info(
            "Tip: Add `APP_BASE_URL` in Streamlit Secrets (so links are clickable). For now, copy this query string and paste after your app URL:"
        )
        st.code(share_link, language="text")

    # IG share card download
    st.subheader("📸 Download IG Share Card")

    # Decide which base image to use for card:
    if apod_ok and apod_img is not None:
        base_for_card = apod_img
        footer = "Powered by NASA APOD + Astropy"
    else:
        base_for_card = generate_starfield_ig(seed_text=f"{birth_date.isoformat()}|{birth_city}|fallback")
        footer = "APOD degraded — fallback visual. Planets still accurate."

    card_buf = create_share_card(
        base_img=base_for_card,
        birth_date=birth_date,
        city=birth_city,
        mode=mode_key,
        visible_planets=visible_planets,
        apod_title=apod_title if apod_ok else "",
        footer_note=footer,
    )

    st.download_button(
        "Download PNG",
        data=card_buf,
        file_name="nasa_birthday_sky.png",
        mime="image/png",
        use_container_width=True,
    )

    # Optional: auto-set query params so refresh keeps current result
    st.experimental_set_query_params(
        date=birth_date.isoformat(),
        city=birth_city,
        mode=mode_key,
        time=exact_time_str,
        autorun="1",
    )
else:
    st.caption("Enter your birthday + city, then press **Run my birthday** 🚀")
