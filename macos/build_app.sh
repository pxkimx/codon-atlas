#!/bin/bash
# Builds "Codon Atlas.app", a double-clickable app with an icon. The program code is
# copied INTO the app (Contents/Resources/app): macOS does not let an app read scripts from
# ~/Documents, ~/Desktop or ~/Downloads without a privacy grant, so an app that pointed back at the
# repository folder would silently fail to start. Rebuild the app after changing the code.
#   macos/build_app.sh [destination-folder]      (default: the repository folder)
set -eu
HERE="$(cd "$(dirname "$0")/.." && pwd)"
DEST="${1:-$HERE}"
APP="$DEST/Codon Atlas.app"
VERSION="$(sed -n 's/^__version__ = "\(.*\)"/\1/p' "$HERE/codon_atlas/__init__.py")"
rm -rf "$APP"; mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
mkdir -p "$APP/Contents/Resources/app"
cp "$HERE/launch.sh" "$HERE/requirements.txt" "$APP/Contents/Resources/app/"
rsync -a --exclude '__pycache__' --exclude '*.pyc' "$HERE/codon_atlas" "$APP/Contents/Resources/app/"
cat > "$APP/Contents/MacOS/CodonAtlas" <<'SH'
#!/bin/bash
export MODE=gui
exec /bin/bash "$(cd "$(dirname "$0")/../Resources/app" && pwd)/launch.sh"
SH
chmod +x "$APP/Contents/MacOS/CodonAtlas"
cat > "$APP/Contents/Info.plist" <<PL
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
  <key>CFBundleName</key><string>Codon Atlas</string>
  <key>CFBundleDisplayName</key><string>Codon Atlas</string>
  <key>CFBundleIdentifier</key><string>io.github.codon-atlas</string>
  <key>CFBundleVersion</key><string>$VERSION</string>
  <key>CFBundleShortVersionString</key><string>$VERSION</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleExecutable</key><string>CodonAtlas</string>
  <key>CFBundleIconFile</key><string>AppIcon</string>
  <key>LSMinimumSystemVersion</key><string>12.0</string>
  <key>NSHighResolutionCapable</key><true/>
</dict></plist>
PL
PY="$(command -v python3)"
"$PY" "$HERE/macos/make_icon.py" "$APP/Contents/Resources/AppIcon.icns"
touch "$APP"
echo "Built: $APP"
