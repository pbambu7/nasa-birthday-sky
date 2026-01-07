import io
import math
import time
import random
from datetime import date, datetime, time as dtime

import streamlit as st
import requests
from PIL import Image, ImageDraw, ImageFont

import pytz
from timezonefinder import TimezoneFinder

from astropy.time import Time
from astropy.coordinates import (
    EarthLocation,
    AltAz,
    get_body,
    SkyCoord,
    BarycentricTrueEcliptic,
)
from astropy import units as u

from geopy.geocoders import Nominatim


# =========================
# CONFIG / PAGE
# =========================
st.set_page_config(page_title="NASA Birthday Sky", layout="wide")

st.markdown(
    """
<style>
.block-container { padding-top: 2rem; padding-bottom: 3rem; max-width: 900px; }
.stButton>button { width:100%; padding:14px; border-radius:14px; font-size:16px; }
small.muted { color: rgba(255,255,255,0.65); }
hr { margin: 1.25rem 0; }
</style>
""",
    unsafe_allow_html=True,
)

NASA_API_KEY = st.secrets.get("NASA_API_KEY")
APP_BASE_URL = st.secrets.get("APP_BASE_URL", "").strip()

if not NASA_API_KEY:
    st.error('NASA API key missing. Add it in Streamlit → App settings → Secrets:\n\nNASA_API_KEY = "YOUR_KEY"')
    st.stop()

st.title("🌌 NASA Birthday Sky")
st.caption("Personal astronomy using NASA open data + real sky calculations (timezone-correct).")


# =========================
# LIGHTWEIGHT "REAL STAR CATALOG" (BRIGHT STAR SUBSET)
# =========================
# NOTE: This is a compact embedded subset (bright stars) with real RA/Dec coordinates.
# It’s “real catalog data” in the sense that these are actual star coordinates,
# but it’s intentionally small so it runs fast + reliably on Streamlit Cloud.
BRIGHT_STARS = {
    # Orion region
    "Betelgeuse": (88.7929, 7.4071),
    "Rigel": (78.6345, -8.2016),
    "Bellatrix": (81.2828, 6.3497),
    "Saiph": (86.9391, -9.6696),
    "Alnitak": (85.1897, -1.9426),
    "Alnilam": (84.0534, -1.2019),
    "Mintaka": (83.0017, -0.2991),

    # Ursa Major (Big Dipper)
    "Dubhe": (165.9322, 61.7510),
    "Merak": (165.4603, 56.3824),
    "Phecda": (178.4577, 53.6948),
    "Megrez": (183.8565, 57.0326),
    "Alioth": (193.5073, 55.9598),
    "Mizar": (200.9814, 54.9254),
    "Alkaid": (206.8856, 49.3133),

    # Leo (a few)
    "Regulus": (152.0929, 11.9672),
    "Denebola": (177.2649, 14.5721),
    "Algieba": (154.9931, 19.8415),

    # Summer triangle
    "Vega": (279.2347, 38.7837),
    "Deneb": (310.3579, 45.2803),
    "Altair": (297.6958, 8.8683),

    # Others
    "Sirius": (101.2875, -16.7161),
    "Canopus": (95.9879, -52.6957),
    "Arcturus": (213.9153, 19.1824),
    "Spica": (201.2983, -11.1613),
    "Antares": (247.3519, -26.4320),
    "Capella": (79.1723, 45.9979),
}

# Constellation line segments defined by star names in BRIGHT_STARS
CONSTELLATION_LINES = {
    "Orion": [
        ("Betelgeuse", "Bellatrix"),
        ("Bellatrix", "Mintaka"),
        ("Mintaka", "Alnilam"),
        ("Alnilam", "Alnitak"),
        ("Alnitak", "Saiph"),
        ("Saiph", "Rigel"),
        ("Rigel", "Betelgeuse"),
    ],
    "Ursa Major": [
        ("Dubhe", "Merak"),
        ("Merak", "Phecda"),
        ("Phecda", "Megrez"),
        ("Megrez", "Alioth"),
        ("Alioth", "Mizar"),
        ("Mizar", "Alkaid"),
    ],
    "Leo": [
        ("Regulus", "Algieba"),
        ("Algieba", "Denebola"),
    ],
    "Summer Triangle": [
        ("Vega", "Deneb"),
        ("Deneb", "Altair"),
        ("Altair", "Vega"),
    ],
}


# =========================
# QUERY PARAMS (SHARE LINKS)
# =========================
def get_query_params():
    # Streamlit has moved APIs across versions; support both.
    try:
        return dict(st.query_params)
    except Exception:
        try:
            qp = st.experimental_get_query_params()
            return {k: v[0] if isinstance(v, list) and v else v for k, v in qp.items()}
        except Exception:
            return {}

qp = get_query_params()

def qp_get(key, default=""):
    v = qp.get(key, default)
    if isinstance(v, list):
        return v[0] if v else default
    return v if v is not None else default

def parse_iso_date(s, fallback):
    try:
        return date.fromisoformat(str(s))
    except Exception:
        return fallback

def parse_hhmm(s, fallback_time=dtime(0, 0)):
    try:
        parts = str(s).split(":")
        hh = max(0, min(23, int(parts[0])))
        mm = max(0, min(59, int(parts[1]) if len(parts) > 1 else 0))
        return dtime(hh, mm)
    except Exception:
        return fallback_time


default_birth_date = parse_iso_date(qp_get("date", "2000-01-01"), date(2000, 1, 1))
default_city = qp_get("city", "")
default_mode = qp_get("mode", "midnight").lower()  # midnight/noon/exact
default_time = parse_hhmm(qp_get("time", "00:00"), dtime(0, 0))
default_autorun = str(qp_get("autorun", "0")) == "1"


# =========================
# UI INPUTS (MOBILE-FRIENDLY)
# =========================
top1, top2 = st.columns(2)
with top1:
    birth_date = st.date_input("Your birth date", value=default_birth_date)
with top2:
    birth_city = st.text_input("Birth city", value=default_city, placeholder="Salt Lake City, UT")

mode_labels = {
    "midnight": "Midnight (00:00)",
    "noon": "Noon (12:00)",
    "exact": "Exact time",
}
mode_options = ["Midnight (00:00)", "Noon (12:00)", "Exact time"]
default_mode_label = mode_labels.get(default_mode, "Midnight (00:00)")

moment = st.selectbox("Moment to calculate (local time)", mode_options, index=mode_options.index(default_mode_label))

birth_time = dtime(0, 0)
if moment.startswith("Midnight"):
    birth_time = dtime(0, 0)
elif moment.startswith("Noon"):
    birth_time = dtime(12, 0)
else:
    birth_time = st.time_input("Exact birth time (optional)", value=default_time)

educator_mode = st.toggle("🧪 Educator Mode (explain the science)", value=False)

st.divider()


# =========================
# CACHE HELPERS
# =========================
@st.cache_data(show_spinner=False, ttl=60 * 60 * 24)
def geocode_city(city: str):
    city = city.strip()
    if not city:
        return None
    try:
        geolocator = Nominatim(user_agent="nasa-birthday-sky", timeout=6)
        loc = geolocator.geocode(city)
        if not loc:
            return None
        return float(loc.latitude), float(loc.longitude), str(loc.address)
    except Exception:
        return None


def local_dt_to_utc(birth_date: date, birth_time: dtime, lat: float, lon: float):
    tf = TimezoneFinder()
    tz_name = tf.timezone_at(lat=lat, lng=lon)
    if not tz_name:
        return None, None

    tz = pytz.timezone(tz_name)
    local_dt = tz.localize(datetime.combine(birth_date, birth_time))
    utc_dt = local_dt.astimezone(pytz.utc)
    return tz_name, utc_dt


@st.cache_data(show_spinner=False, ttl=60 * 60)
def fetch_apod(date_iso: str, api_key: str):
    url = "https://api.nasa.gov/planetary/apod"
    params = {"date": date_iso, "api_key": api_key}

    # Robust retries for degraded APOD
    for attempt in range(4):
        try:
            r = requests.get(url, params=params, timeout=(4, 10))
            r.raise_for_status()
            return r.json()
        except requests.exceptions.Timeout:
            time.sleep(0.6 * (2 ** attempt))
        except Exception:
            return None
    return None


# =========================
# ASTRONOMY: PLANETS + MOON
# =========================
def compute_visible_planets(location: EarthLocation, obs_time: Time):
    frame = AltAz(obstime=obs_time, location=location)
    planet_names = ["venus", "mars", "jupiter"]
    results = []

    for p in planet_names:
        try:
            body = get_body(p, obs_time).transform_to(frame)
            alt = body.alt.to(u.deg).value
            az = body.az.to(u.deg).value
            if alt > 0:
                results.append({"name": p.capitalize(), "alt": float(alt), "az": float(az)})
        except Exception:
            continue

    results.sort(key=lambda x: x["alt"], reverse=True)
    return results


def compute_moon_phase(obs_time: Time):
    """
    Returns:
    - illumination fraction (0..1)
    - phase name (string)
    - waxing bool
    Uses Sun–Moon elongation in the sky:
      k = (1 - cos(D)) / 2
    """
    moon = get_body("moon", obs_time)
    sun = get_body("sun", obs_time)

    # Angular separation (elongation) as seen from Earth
    D = moon.separation(sun).to(u.rad).value  # radians
    illum = (1 - math.cos(D)) / 2  # 0 new -> 1 full

    # Waxing/waning: compare ecliptic longitudes
    moon_ecl = moon.transform_to(BarycentricTrueEcliptic(obstime=obs_time))
    sun_ecl = sun.transform_to(BarycentricTrueEcliptic(obstime=obs_time))
    dlon = (moon_ecl.lon - sun_ecl.lon).wrap_at(360 * u.deg).to(u.deg).value
    waxing = dlon > 0

    # Phase name
    pct = illum * 100
    if pct < 2:
        name = "New Moon"
    elif pct < 48:
        name = "Waxing Crescent" if waxing else "Waning Crescent"
    elif pct < 52:
        name = "First Quarter" if waxing else "Last Quarter"
    elif pct < 98:
        name = "Waxing Gibbous" if waxing else "Waning Gibbous"
    else:
        name = "Full Moon"

    return float(illum), name, waxing


# =========================
# VISUALS: SKY + CONSTELLATIONS + PLANETS
# =========================
def altaz_to_xy(alt_deg, az_deg, w, h):
    """
    Simple projection:
    - x based on azimuth 0..360
    - y based on altitude 0..90 (top)
    """
    x = int((az_deg / 360.0) * w)
    # keep horizon near bottom; clamp alt for view
    alt_clamped = max(0.0, min(90.0, float(alt_deg)))
    y = int(h - (alt_clamped / 90.0) * (h * 0.85) - (h * 0.05))
    return x, y


def draw_glow(draw, x, y, base_radius=4, glow_radius=14, color=(255, 215, 0)):
    # glow
    for r in range(glow_radius, base_radius, -2):
        alpha = int(18 * (r / glow_radius))
        # simulate alpha by blending toward black (simple hack)
        c = tuple(int(color[i] * (alpha / 255)) for i in range(3))
        draw.ellipse((x - r, y - r, x + r, y + r), fill=c)
    # core
    draw.ellipse((x - base_radius, y - base_radius, x + base_radius, y + base_radius), fill=color)


def project_star_to_altaz(star_ra_deg, star_dec_deg, obs_time: Time, location: EarthLocation):
    sc = SkyCoord(ra=star_ra_deg * u.deg, dec=star_dec_deg * u.deg, frame="icrs")
    altaz = sc.transform_to(AltAz(obstime=obs_time, location=location))
    return altaz.alt.to(u.deg).value, altaz.az.to(u.deg).value


def generate_sky_image(
    obs_time: Time,
    location: EarthLocation,
    visible_planets,
    show_constellations=True,
    caption_text="Simulated night sky based on your birth location",
):
    W, H = 1080, 620
    img = Image.new("RGB", (W, H), (0, 0, 0))
    draw = ImageDraw.Draw(img)

    # starfield background (random but stable-ish per day)
    seed = int(obs_time.jd * 10) % 10_000_000
    rng = random.Random(seed)

    # stars
    for _ in range(1400):
        x = rng.randint(0, W - 1)
        y = rng.randint(0, H - 1)
        b = int(80 + (rng.random() ** 0.4) * 175)
        r = 1 if rng.random() < 0.95 else 2
        draw.ellipse((x - r, y - r, x + r, y + r), fill=(b, b, b))

    # constellations
    if show_constellations:
        line_color = (90, 140, 255)  # subtle
        label_color = (150, 180, 255)

        for cname, segs in CONSTELLATION_LINES.items():
            # draw each segment if both stars are above horizon
            for a, b in segs:
                if a not in BRIGHT_STARS or b not in BRIGHT_STARS:
                    continue

                a_ra, a_dec = BRIGHT_STARS[a]
                b_ra, b_dec = BRIGHT_STARS[b]

                a_alt, a_az = project_star_to_altaz(a_ra, a_dec, obs_time, location)
                b_alt, b_az = project_star_to_altaz(b_ra, b_dec, obs_time, location)

                if a_alt <= 0 and b_alt <= 0:
                    continue

                ax, ay = altaz_to_xy(a_alt, a_az, W, H)
                bx, by = altaz_to_xy(b_alt, b_az, W, H)
                draw.line((ax, ay, bx, by), fill=line_color, width=1)

            # label: pick first star in first segment if visible
            first_seg = segs[0] if segs else None
            if first_seg:
                sname = first_seg[0]
                if sname in BRIGHT_STARS:
                    ra, dec = BRIGHT_STARS[sname]
                    alt, az = project_star_to_altaz(ra, dec, obs_time, location)
                    if alt > 0:
                        x, y = altaz_to_xy(alt, az, W, H)
                        draw.text((x + 8, y - 12), cname, fill=label_color)

    # planets markers
    for p in visible_planets:
        x, y = altaz_to_xy(p["alt"], p["az"], W, H)
        draw_glow(draw, x, y, base_radius=4, glow_radius=16, color=(255, 215, 0))
        draw.text((x + 10, y - 10), p["name"], fill=(255, 255, 255))

    # bottom overlay caption (Apple/NASA feel)
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    o = ImageDraw.Draw(overlay)
    o.rectangle((0, H - 80, W, H), fill=(0, 0, 0, 160))
    o.text((24, H - 58), caption_text, fill=(255, 255, 255))
    img = Image.alpha_composite(img.convert("RGBA"), overlay).convert("RGB")

    return img


def create_share_card_ig(base_img: Image.Image, birth_date: date, city: str, moment_label: str, tz_name: str,
                         planets_text: str, moon_text: str):
    # 1080 x 1080
    img = base_img.convert("RGB")

    # make square (center crop or pad)
    if img.width != img.height:
        # center-crop to square
        side = min(img.width, img.height)
        left = (img.width - side) // 2
        top = (img.height - side) // 2
        img = img.crop((left, top, left + side, top + side))
    img = img.resize((1080, 1080))

    draw = ImageDraw.Draw(img)
    draw.rectangle((0, 760, 1080, 1080), fill=(0, 0, 0))

    # fonts (fallback safe)
    try:
        font_big = ImageFont.truetype("DejaVuSans.ttf", 44)
        font_mid = ImageFont.truetype("DejaVuSans.ttf", 30)
        font_small = ImageFont.truetype("DejaVuSans.ttf", 24)
    except Exception:
        font_big = font_mid = font_small = None

    draw.text((40, 790), "NASA Birthday Sky", fill="white", font=font_big)
    draw.text((40, 850), f"{birth_date.isoformat()} • {city}", fill="white", font=font_mid)
    draw.text((40, 895), f"{moment_label} • {tz_name}", fill=(200, 200, 200), font=font_small)
    draw.text((40, 935), f"Planets: {planets_text}", fill=(200, 200, 200), font=font_small)
    draw.text((40, 970), f"Moon: {moon_text}", fill=(200, 200, 200), font=font_small)
    draw.text((40, 1030), "Powered by NASA Open APIs • Astropy", fill=(150, 150, 150), font=font_small)

    buf = io.BytesIO()
    img.save(buf, format="PNG")
    buf.seek(0)
    return buf


def build_share_url(birth_date: date, city: str, mode_key: str, time_hhmm: str):
    q_city = city.replace(" ", "%20")
    qs = f"?date={birth_date.isoformat()}&city={q_city}&mode={mode_key}&time={time_hhmm}&autorun=1"
    if APP_BASE_URL:
        return APP_BASE_URL.rstrip("/") + "/" + qs
    return qs


# =========================
# RUN ACTION
# =========================
run_clicked = st.button("🚀 Run my birthday")
should_autorun = default_autorun
should_run = run_clicked or should_autorun

if should_run:
    if not birth_city.strip():
        st.warning("Please enter a birth city (example: “Kinshasa, Congo” or “Salt Lake City, UT”).")
        st.stop()

    with st.spinner("Finding your city coordinates..."):
        geo = geocode_city(birth_city)

    if not geo:
        st.error("Couldn’t find that city. Try adding country/state (e.g., “Kinshasa, Congo”).")
        st.stop()

    lat, lon, resolved_address = geo
    location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg)

    # mode key
    if moment.startswith("Midnight"):
        mode_key = "midnight"
        moment_label = "Midnight"
        birth_time_final = dtime(0, 0)
    elif moment.startswith("Noon"):
        mode_key = "noon"
        moment_label = "Noon"
        birth_time_final = dtime(12, 0)
    else:
        mode_key = "exact"
        moment_label = "Exact time"
        birth_time_final = birth_time

    tz_name, utc_dt = local_dt_to_utc(birth_date, birth_time_final, lat, lon)
    if not tz_name or not utc_dt:
        st.error("Could not determine timezone for this location.")
        st.stop()

    obs_time = Time(utc_dt)

    # compute planets + moon
    with st.spinner("Calculating planets + moon (timezone-correct)..."):
        visible_planets = compute_visible_planets(location, obs_time)
        moon_illum, moon_phase_name, waxing = compute_moon_phase(obs_time)

    moon_pct = int(round(moon_illum * 100))
    moon_text = f"{moon_phase_name} • {moon_pct}% lit"

    # nice animated planet sweep
    st.subheader("🪐 Planets visible above the horizon")
    if visible_planets:
        ph = st.empty()
        lines = []
        for p in visible_planets:
            lines.append(f"**{p['name']}** — altitude **{p['alt']:.1f}°**, azimuth **{p['az']:.1f}°**")
            ph.markdown("\n\n".join(lines))
            time.sleep(0.35)
    else:
        st.write("No major planets were above the horizon at that moment (Mars/Venus/Jupiter list).")

    st.subheader("🌕 Moon")
    st.write(f"**{moon_phase_name}** — **{moon_pct}%** illuminated")

    # Educator mode explanation
    if educator_mode:
        with st.expander("🧪 How the calculations work (Educator Mode)", expanded=True):
            st.markdown(
                f"""
**1) Location → Coordinates**  
We convert your city into latitude/longitude using geocoding.

**2) Coordinates → Timezone → UTC**  
Your input time is local to the birth city.  
We determine the timezone (**{tz_name}**) and convert local time to **UTC** (the standard for astronomy).

**3) Planet visibility = altitude > 0°**  
Altitude tells how high an object is above the horizon:
- **Altitude > 0°** → visible above horizon  
- **Altitude ≤ 0°** → below horizon

**4) Azimuth is direction**  
Azimuth is a compass direction in degrees:
- 0° = North, 90° = East, 180° = South, 270° = West

**5) Moon illumination**  
We compute Sun–Moon separation and estimate illumination:
- New Moon ≈ 0%
- Full Moon ≈ 100%
"""
            )

    # NASA APOD (optional; may be degraded)
    with st.spinner("Fetching NASA APOD (may be degraded during outages)..."):
        apod = fetch_apod(birth_date.isoformat(), NASA_API_KEY)

    apod_img = None
    apod_title = ""
    apod_expl = ""
    apod_ok = False

    if isinstance(apod, dict):
        apod_title = str(apod.get("title", "")).strip()
        apod_expl = str(apod.get("explanation", "")).strip()
        apod_url = str(apod.get("url", "")).strip()
        media_type = str(apod.get("media_type", "image")).lower()

        if media_type == "image" and apod_url:
            try:
                rr = requests.get(apod_url, timeout=(4, 10))
                rr.raise_for_status()
                apod_img = Image.open(io.BytesIO(rr.content)).convert("RGB")
                apod_ok = True
            except Exception:
                apod_ok = False

    st.divider()

    # Visual section
    if apod_ok and apod_img is not None:
        st.success("✅ NASA APOD loaded successfully.")
        st.subheader(apod_title or "NASA Astronomy Picture of the Day")
        st.image(apod_img, use_container_width=True)
        if apod_expl:
            with st.expander("APOD explanation"):
                st.write(apod_expl)
        visual_for_card = apod_img
    else:
        st.warning("🚨 NASA APOD is temporarily unavailable. Showing a simulated sky visualization (planets remain accurate).")
        sky_img = generate_sky_image(
            obs_time=obs_time,
            location=location,
            visible_planets=visible_planets,
            show_constellations=True,
            caption_text="Simulated night sky based on your birth location",
        )
        st.image(sky_img, use_container_width=True, caption="Simulated night sky based on your birth location")
        visual_for_card = sky_img

    # Share link + download card
    planets_text = ", ".join([p["name"] for p in visible_planets]) if visible_planets else "None"
    time_hhmm = f"{birth_time_final.hour:02d}:{birth_time_final.minute:02d}"

    st.subheader("🔗 Share this result")
    share_url = build_share_url(birth_date, birth_city, mode_key, time_hhmm)
    st.code(share_url, language="text")
    if not APP_BASE_URL:
        st.caption("Tip: add APP_BASE_URL in Streamlit Secrets for a full clickable link.")

    st.subheader("📸 Download IG Share Card")
    card_buf = create_share_card_ig(
        base_img=visual_for_card,
        birth_date=birth_date,
        city=birth_city,
        moment_label=moment_label,
        tz_name=tz_name,
        planets_text=planets_text,
        moon_text=moon_text,
    )
    st.download_button(
        "Download PNG",
        data=card_buf,
        file_name="nasa_birthday_sky.png",
        mime="image/png",
        use_container_width=True,
    )

    st.caption(f"📍 Location used: {resolved_address} (lat {lat:.3f}, lon {lon:.3f}) • Time used (UTC): {utc_dt.isoformat()}")

    # Persist params so refresh keeps result
    try:
        st.query_params.update({
            "date": birth_date.isoformat(),
            "city": birth_city,
            "mode": mode_key,
            "time": time_hhmm,
            "autorun": "1",
        })
    except Exception:
        st.experimental_set_query_params(
            date=birth_date.isoformat(),
            city=birth_city,
            mode=mode_key,
            time=time_hhmm,
            autorun="1",
        )

else:
    st.caption("Enter your birthday + city, then press **Run my birthday** 🚀")
