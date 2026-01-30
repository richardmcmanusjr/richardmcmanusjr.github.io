# Image Loading & Caching Optimization Guide

## What We've Implemented

Your website now has three layers of image optimization to eliminate lag and improve performance:

### 1. **Service Worker (`sw.js`)**
- **What it does**: Caches all images, CSS, JavaScript, and HTML pages on the user's device
- **Benefits**: 
  - First visit: Normal loading
  - Subsequent visits: Lightning-fast load times (images served from cache)
  - Works offline or with slow connections
  - Automatically updates cache when new files are deployed
- **How it works**: Every time a user visits, the browser stores a copy locally. Next time they visit, it loads the cached version first.

### 2. **Image Preloader (`image-preloader.js`)**
- **What it does**: Preloads all images in batches in the background after the page loads
- **Benefits**:
  - Carousel images are ready instantly when users interact with them
  - Eliminates the "spinning wheel" when switching to different activity/team images
  - Doesn't slow down initial page load (runs in background)
- **How it works**: Loads images in groups of 5 with small delays between batches to avoid overloading the browser

### 3. **HTML Link Preload Tags**
- **What it does**: Tells the browser which images are critical and should load first
- **Benefits**: Critical images (team photos, research icons) start loading as soon as the HTML is parsed
- **Preloaded images per page**:
  - **index.html**: Team member photo, research icons
  - **team.html**: All main team member photos
  - **research.html**: Research category icons
  - **activities.html**: First image from each activity carousel

## Performance Impact

### Before Optimization
- Images loaded on-demand (when carousel scrolled or page section reached)
- Lag/delay when interacting with carousels
- No offline access

### After Optimization
- ✅ Carousel images preloaded and ready instantly
- ✅ All images cached locally after first visit
- ✅ Extremely fast page loads on repeat visits
- ✅ Offline support (cached content loads if network unavailable)
- ✅ Reduced server bandwidth (cached images don't re-download)

## Browser Support

- **Service Workers**: Supported in all modern browsers (Chrome, Firefox, Safari 11.1+, Edge)
- **Image Preloading**: Supported everywhere
- **Graceful Degradation**: If a browser doesn't support service workers, the site still works normally

## Deploying to GitHub Pages

Since your site is on GitHub Pages (`richardmcmanusjr.github.io`), the service worker will:
1. Automatically activate on first visit
2. Persist across page refreshes
3. Auto-update whenever you push new files to GitHub

## Monitoring Cache

To verify the cache is working:
1. Open your site in Chrome/Firefox
2. Open DevTools (F12 or Cmd+Option+I)
3. Go to **Application** tab → **Service Workers**
4. You should see "Service Worker registered successfully" in the console
5. Go to **Application** → **Cache Storage** → **bad-lab-v1** to see cached files

## Adding New Images

If you add new images to your site:

1. Add them to the `imagesToPreload` array in `image-preloader.js`
2. Add them to the `URLS_TO_CACHE` array in `sw.js`
3. Push changes to GitHub

The cache will automatically update on the next visit.

## Troubleshooting

**Images not loading?**
- Clear the cache: DevTools → Application → Cache Storage → bad-lab-v1 → Delete
- Hard refresh the page (Cmd+Shift+R or Ctrl+Shift+F5)

**Service Worker not registering?**
- Check console for errors (F12)
- Ensure `sw.js` is in the root directory
- Service workers only work on HTTPS (your GitHub Pages site uses HTTPS ✓)

## Additional Optimization Tips

To further improve performance:

1. **Image Compression**: Consider compressing JPG/PNG files (use tools like TinyPNG or ImageOptim)
2. **WebP Format**: Convert images to WebP for better compression (85% smaller files)
3. **Lazy Loading**: For very large image galleries, add `loading="lazy"` to `<img>` tags
4. **CDN**: For very large sites, consider hosting images on a CDN

## Technical Details

### Service Worker Scope
The service worker intercepts all network requests from your site and:
- Returns cached version if available
- Falls back to network if not cached
- Caches new responses for future visits
- Handles offline scenarios gracefully

### Cache Strategy: Cache First
```
User Request → Check Cache → If Found: Return Cache
                          → If Not Found: Fetch from Network → Cache & Return
```

This is optimal for images and static assets that change infrequently.

---

**Questions?** Check the browser console (F12) for any error messages, or review the scripts in `sw.js` and `image-preloader.js` for details.
