// Video loading state management
document.addEventListener('DOMContentLoaded', () => {
    const videos = document.querySelectorAll('.hero-video');
    const loadingOverlay = document.getElementById('heroLoadingOverlay');
    let videosReady = 0;
    const totalVideos = videos.length;
    
    function checkAllVideosReady() {
        if (videosReady === totalVideos && loadingOverlay) {
            // All videos are ready, hide the loading overlay
            loadingOverlay.classList.add('loaded');
        }
    }
    
    videos.forEach((video) => {
        // Mark video as ready when it can play
        video.addEventListener('canplay', () => {
            videosReady++;
            checkAllVideosReady();
        }, { once: true });
        
        // Also check if video is already playing
        if (video.readyState >= 2) {
            videosReady++;
            checkAllVideosReady();
        }
    });
    
    // Fallback: hide loading after max time to ensure page remains functional
    setTimeout(() => {
        if (loadingOverlay && !loadingOverlay.classList.contains('loaded')) {
            loadingOverlay.classList.add('loaded');
        }
    }, 5000);
});

// Mobile coach image blur/focus on scroll
document.addEventListener('DOMContentLoaded', () => {
    function isMobile() {
        return window.innerWidth <= 768;
    }
    const coachImg = document.querySelector('.coach-image img');
    if (coachImg) {
        function handleMobileCoachImage() {
            if (!isMobile()) {
                coachImg.classList.remove('mobile-blur', 'mobile-focus');
                return;
            }
            const coachSection = document.querySelector('.coach-section');
            const rect = coachImg.getBoundingClientRect();
            const windowHeight = window.innerHeight;
            const imgCenter = rect.top + rect.height / 2;
            const windowCenter = windowHeight / 2;
            const distance = Math.abs(imgCenter - windowCenter);
            
            // Check if coach section is in view at all
            if (rect.bottom < 0 || rect.top > windowHeight) {
                coachImg.classList.remove('mobile-focus');
                coachImg.classList.add('mobile-blur');
                return;
            }
            
            // Focus if image center is close to window center (50% of window height threshold for mobile)
            // This prevents blur when scrolling past the section
            if (distance < windowHeight * 0.25) {
                coachImg.classList.add('mobile-focus');
                coachImg.classList.remove('mobile-blur');
            } else {
                coachImg.classList.add('mobile-blur');
                coachImg.classList.remove('mobile-focus');
            }
        }
        // Use passive listeners for better scroll performance
        window.addEventListener('scroll', handleMobileCoachImage, { passive: true });
        window.addEventListener('resize', handleMobileCoachImage);
        handleMobileCoachImage();
    }
});
// GOATCoach popup animation for 'history of speed' hover and click
document.addEventListener('DOMContentLoaded', () => {
    const history = document.querySelector('.history-of-speed');
    const popup = document.querySelector('.goatcoach-popup');
    if (history && popup) {
        let hideTimeout;
        
        const showPopup = () => {
            popup.style.opacity = '1';
            popup.style.pointerEvents = 'auto';
            popup.style.transform = 'scale(1.1) rotate(0deg) translateY(-10px)';
        };
        
        const hidePopup = () => {
            popup.style.opacity = '';
            popup.style.pointerEvents = '';
            popup.style.transform = '';
        };
        
        const scheduleHide = () => {
            clearTimeout(hideTimeout);
            hideTimeout = setTimeout(hidePopup, 3000);
        };
        
        // Desktop hover events
        history.addEventListener('mouseenter', showPopup);
        history.addEventListener('mouseleave', scheduleHide);
        
        // Focus/blur for keyboard navigation
        history.addEventListener('focus', showPopup);
        history.addEventListener('blur', hidePopup);
        
        // Mobile click/touch support - show for 3 seconds
        history.addEventListener('click', (e) => {
            e.stopPropagation();
            showPopup();
            scheduleHide();
        });
    }
});

// Charles popup animation for 'Charles the turkey' hover and click
document.addEventListener('DOMContentLoaded', () => {
    const charles = document.querySelector('.charles-turkey');
    const popup = document.querySelector('.charles-popup');
    if (charles && popup) {
        let hideTimeout;
        
        const showPopup = () => {
            popup.style.opacity = '1';
            popup.style.pointerEvents = 'auto';
            popup.style.transform = 'scale(1.1) rotate(0deg) translateY(-10px)';
        };
        
        const hidePopup = () => {
            popup.style.opacity = '';
            popup.style.pointerEvents = '';
            popup.style.transform = '';
        };
        
        const scheduleHide = () => {
            clearTimeout(hideTimeout);
            hideTimeout = setTimeout(hidePopup, 3000);
        };
        
        // Desktop hover events
        charles.addEventListener('mouseenter', showPopup);
        charles.addEventListener('mouseleave', hidePopup);
        
        // Focus/blur for keyboard navigation
        charles.addEventListener('focus', showPopup);
        charles.addEventListener('blur', hidePopup);
        
        // Mobile click/touch support - show for 3 seconds
        charles.addEventListener('click', (e) => {
            e.stopPropagation();
            showPopup();
            scheduleHide();
        });
    }
});
// Apple-style scroll animations
document.addEventListener('DOMContentLoaded', () => {
    // Smooth scroll observer for elements
    const observerOptions = {
        threshold: 0.15,
        rootMargin: '0px 0px -100px 0px'
    };

    const observer = new IntersectionObserver((entries) => {
        entries.forEach((entry, index) => {
            if (entry.isIntersecting) {
                // Add staggered delay for multiple elements
                setTimeout(() => {
                    entry.target.classList.add('visible');
                }, index * 100);
            }
        });
    }, observerOptions);

    // Observe all animatable elements
    const animatableElements = document.querySelectorAll(
        '.stat-card, .record-card, .philosophy-item, .records h2, .philosophy-content h2, .featured-image, .featured-text'
    );
    
    animatableElements.forEach(el => {
        observer.observe(el);
    });

    // After initial load animations finish, clear them so scroll transforms apply
    let animationTimeoutId;
    const clearAnimationTimeout = () => {
        clearTimeout(animationTimeoutId);
    };
    
    animationTimeoutId = setTimeout(() => {
        const hero = document.querySelector('.hero-content');
        const title = document.querySelector('.hero h1');
        const subtitle = document.querySelector('.hero-subtitle');
        [hero, title, subtitle].forEach(el => {
            if (!el) return;
            el.style.animation = 'none';
            el.style.opacity = '1';
            el.style.transform = 'none';
            el.style.filter = 'none';
        });
    }, 2000);

    // Parallax effect removed - hero stays fixed
    // Add scroll effects to hero text and coach section
    let ticking = false;
    let heroAnimationDisabled = false;
    window.addEventListener('scroll', () => {
        clearAnimationTimeout();
        
        if (!ticking) {
            window.requestAnimationFrame(() => {
                const scrolled = window.pageYOffset;
                const hero = document.querySelector('.hero-content');
                const backLink = document.querySelector('.back-link');
                const title = document.querySelector('.hero h1');
                const subtitle = document.querySelector('.hero-subtitle');
                const heroOverlay = document.querySelector('.hero-overlay');
                const roadmapSection = document.querySelector('.roadmap');
                
                // Check if timeline section is in view
                let timelineInView = false;
                if (roadmapSection) {
                    const roadmapRect = roadmapSection.getBoundingClientRect();
                    timelineInView = roadmapRect.top < window.innerHeight * 0.5;
                }
                
                if (hero) {
                    if (scrolled >= window.innerHeight * 0.8 || timelineInView) {
                        // Force complete hide when scrolled past threshold or timeline is visible
                        // Only remove animations once when actually scrolling down
                        if (!heroAnimationDisabled) {
                            hero.style.animation = 'none';
                            if (title) title.style.animation = 'none';
                            if (subtitle) subtitle.style.animation = 'none';
                            heroAnimationDisabled = true;
                        }
                        hero.style.opacity = '0';
                        hero.style.pointerEvents = 'none';
                        hero.style.transform = 'none';
                        hero.style.filter = 'none';
                    } else if (scrolled < window.innerHeight * 1.5) {
                        // Back near the top - restore hero visibility
                        // Re-enable animations when scrolling back near the top
                        if (heroAnimationDisabled && scrolled < window.innerHeight * 0.3) {
                            heroAnimationDisabled = false;
                            if (title) title.style.animation = '';
                            if (subtitle) subtitle.style.animation = '';
                            hero.style.animation = '';
                        }
                        hero.style.pointerEvents = 'auto';
                        const opacity = 1 - (scrolled / (window.innerHeight * 0.8));
                        const scale = 1 - (scrolled / (window.innerHeight * 2.2));
                        const blur = (scrolled / 50);
                        const translateY = -(scrolled * 0.08);
                        
                        hero.style.opacity = Math.max(0, opacity);
                        hero.style.transform = `scale(${Math.max(0.6, scale)}) translateY(${translateY}px)`;
                        hero.style.filter = `blur(${Math.min(blur, 10)}px)`;
                    }
                }
                
                // Fade back-link along with hero
                if (backLink) {
                    const opacity = 1 - (scrolled / (window.innerHeight * 0.8));
                    backLink.style.opacity = Math.max(0, opacity);
                }
                
                if (heroOverlay) {
                    if (scrolled >= window.innerHeight || timelineInView) {
                        heroOverlay.style.opacity = '0';
                    } else {
                        const overlayOpacity = 1 - (scrolled / (window.innerHeight * 1.2));
                        heroOverlay.style.opacity = Math.max(0, overlayOpacity);
                    }
                }

                // Coach section - desktop/mobile responsive handling
                const coachSection = document.querySelector('.coach-section');
                const coachCard = document.querySelector('.coach-card');
                const coachImage = document.querySelector('.coach-image');
                const isMobile = window.innerWidth <= 768;
                
                if (coachSection && isMobile) {
                    // MOBILE: Fixed position, centered, with fade-in/out based on scroll
                    const roadmapSection = document.querySelector('.roadmap');
                    const featuredSection = document.querySelector('.featured-race');
                    
                    let coachScrollStart = window.innerHeight * 0.8;
                    let coachScrollEnd = window.innerHeight * 3.2;
                    
                    if (roadmapSection && featuredSection) {
                        const roadmapRect = roadmapSection.getBoundingClientRect();
                        const featuredRect = featuredSection.getBoundingClientRect();
                        
                        const roadmapBottom = roadmapRect.bottom + window.pageYOffset;
                        const featuredTop = featuredRect.top + window.pageYOffset;
                        
                        coachScrollStart = roadmapBottom - (window.innerHeight * 0.6);
                        coachScrollEnd = featuredTop + (window.innerHeight * 0.6);
                    }
                    
                    const transitionZone = window.innerHeight * 0.4;
                    
                    if (scrolled < coachScrollStart) {
                        coachSection.style.opacity = '0';
                        coachSection.style.pointerEvents = 'none';
                    } else if (scrolled >= coachScrollStart && scrolled < coachScrollStart + transitionZone) {
                        coachSection.style.pointerEvents = 'auto';
                        const fadeInProgress = (scrolled - coachScrollStart) / transitionZone;
                        coachSection.style.opacity = Math.min(1, fadeInProgress);
                    } else if (scrolled >= coachScrollStart + transitionZone && scrolled < coachScrollEnd - transitionZone) {
                        coachSection.style.pointerEvents = 'auto';
                        coachSection.style.opacity = '1';
                    } else if (scrolled >= coachScrollEnd - transitionZone && scrolled < coachScrollEnd) {
                        coachSection.style.pointerEvents = 'auto';
                        const fadeOutProgress = (scrolled - (coachScrollEnd - transitionZone)) / transitionZone;
                        coachSection.style.opacity = Math.max(0, 1 - fadeOutProgress);
                    } else if (scrolled >= coachScrollEnd) {
                        coachSection.style.opacity = '0';
                        coachSection.style.pointerEvents = 'none';
                    }
                } else if (coachSection && !isMobile) {
                    // DESKTOP: Fixed position with scroll-based animations
                    const roadmapSection = document.querySelector('.roadmap');
                    const featuredSection = document.querySelector('.featured-race');
                    
                    let coachScrollStart = window.innerHeight * 0.8;
                    let coachScrollEnd = window.innerHeight * 3.2;
                    
                    if (roadmapSection && featuredSection) {
                        const roadmapRect = roadmapSection.getBoundingClientRect();
                        const featuredRect = featuredSection.getBoundingClientRect();
                        
                        const roadmapBottom = roadmapRect.bottom + window.pageYOffset;
                        const featuredTop = featuredRect.top + window.pageYOffset;
                        
                        coachScrollStart = roadmapBottom - (window.innerHeight * 0.6);
                        coachScrollEnd = featuredTop + (window.innerHeight * 0.6);
                    }
                    
                    const transitionZone = window.innerHeight * 0.4;
                    
                    if (scrolled < coachScrollStart) {
                        coachSection.style.opacity = '0';
                        coachSection.style.pointerEvents = 'none';
                    } else if (scrolled >= coachScrollStart && scrolled < coachScrollStart + transitionZone) {
                        coachSection.style.pointerEvents = 'auto';
                        const fadeInProgress = (scrolled - coachScrollStart) / transitionZone;
                        coachSection.style.opacity = Math.min(1, fadeInProgress);
                    } else if (scrolled >= coachScrollStart + transitionZone && scrolled < coachScrollEnd - transitionZone) {
                        coachSection.style.pointerEvents = 'auto';
                        coachSection.style.opacity = '1';
                    } else if (scrolled >= coachScrollEnd - transitionZone && scrolled < coachScrollEnd) {
                        coachSection.style.pointerEvents = 'auto';
                        const fadeOutProgress = (scrolled - (coachScrollEnd - transitionZone)) / transitionZone;
                        coachSection.style.opacity = Math.max(0, 1 - fadeOutProgress);
                    } else if (scrolled >= coachScrollEnd) {
                        coachSection.style.opacity = '0';
                        coachSection.style.pointerEvents = 'none';
                    }
                    
                    // Desktop: Keep card centered in viewport
                    if (coachCard && coachImage) {
                        const coachRect = coachCard.getBoundingClientRect();
                        const coachCenter = coachRect.top + coachRect.height / 2;
                        const viewportCenter = window.innerHeight / 2;
                        const distance = coachCenter - viewportCenter;
                        
                        const centeringOffset = -distance * 0.5;
                        const imageParallaxOffset = distance * 0.15;
                        coachImage.style.transform = `translateY(${imageParallaxOffset}px) scale(1.05)`;
                        
                        const tiltAmount = Math.sin(scrolled * 0.002) * 3;
                        coachCard.style.transform = `translateY(${centeringOffset}px) perspective(1200px) rotateX(${tiltAmount}deg)`;
                    }
                }
                
                ticking = false;
            });
            ticking = true;
        }
    });

    // Smooth navbar background on scroll
    const navbar = document.querySelector('.navbar');
    let lastScroll = 0;
    
    if (navbar) {
        window.addEventListener('scroll', () => {
            const currentScroll = window.pageYOffset;
            
            if (currentScroll > 100) {
                navbar.style.background = 'rgba(255, 255, 255, 0.95)';
                navbar.style.backdropFilter = 'blur(20px)';
            } else {
                navbar.style.background = 'rgba(255, 255, 255, 0.8)';
                navbar.style.backdropFilter = 'blur(20px)';
            }
            
            lastScroll = currentScroll;
        });
    }

    // Add micro-interactions to cards
    const cards = document.querySelectorAll('.stat-card, .record-card');
    cards.forEach(card => {
        card.addEventListener('mouseenter', function() {
            this.style.transition = 'all 0.4s cubic-bezier(0.25, 0.46, 0.45, 0.94)';
        });
        
        card.addEventListener('mousemove', function(e) {
            const rect = this.getBoundingClientRect();
            const x = e.clientX - rect.left;
            const y = e.clientY - rect.top;
            
            const centerX = rect.width / 2;
            const centerY = rect.height / 2;
            
            const rotateX = (y - centerY) / 20;
            const rotateY = (centerX - x) / 20;
            
            this.style.transform = `translateY(-15px) scale(1.03) perspective(1000px) rotateX(${rotateX}deg) rotateY(${rotateY}deg)`;
        });
        
        card.addEventListener('mouseleave', function() {
            this.style.transform = '';
        });
    });

    // Smooth scroll for navigation links
    document.querySelectorAll('a[href^="#"]').forEach(anchor => {
        anchor.addEventListener('click', function (e) {
            e.preventDefault();
            const target = document.querySelector(this.getAttribute('href'));
            if (target) {
                window.scrollTo({
                    top: target.offsetTop - 100,
                    behavior: 'smooth'
                });
            }
        });
    });

    // Timeline "Now" clock
    const nowEl = document.getElementById('timeline-now');
    const nowLabel = nowEl ? nowEl.querySelector('.now-label') : null;

    function updateNowClock() {
        if (!nowLabel) return;
        nowLabel.textContent = `Today`;
    }

    updateNowClock();
    // Update every 30 seconds to keep it reasonably fresh
    setInterval(updateNowClock, 30000);

    // Position the now marker relative to event dates along the timeline
    const timeline = document.querySelector('.timeline');
    const timelineLine = document.querySelector('.timeline-line');
    const timelineNow = document.getElementById('timeline-now');

    function parseDateAttr(el) {
        const v = el.getAttribute('data-date');
        if (!v) return null;
        const d = new Date(v + 'T00:00:00');
        return isNaN(d) ? null : d;
    }

    function positionNowMarker() {
        if (!timeline || !timelineLine || !timelineNow) return;
        // collect event dates
        const items = Array.from(timeline.querySelectorAll('.timeline-item')).map(item => ({ el: item, date: parseDateAttr(item) }));
        const dates = items.map(i => i.date).filter(Boolean);
        if (dates.length === 0) return;

        const minDate = new Date(Math.min(...dates.map(d => d.getTime())));
        const maxDate = new Date(Math.max(...dates.map(d => d.getTime())));
        const now = new Date();

        const lineRect = timelineLine.getBoundingClientRect();
        const timelineRect = timeline.getBoundingClientRect();

        // compute fraction between min and max
        const span = maxDate.getTime() - minDate.getTime();
        let frac = 0;
        if (span > 0) {
            frac = (now.getTime() - minDate.getTime()) / span;
            frac = Math.max(0, Math.min(1, frac));
        }

        // respect responsive orientation (vertical timeline on small screens)
        const isNarrow = window.matchMedia('(max-width: 1024px)').matches;
        if (!isNarrow) {
            // horizontal: compute left in px relative to timeline container
            // clear minHeight that was set for narrow screens
            timeline.style.minHeight = '';
            
            const lineLeft = lineRect.left - timelineRect.left;
            const leftPx = lineLeft + (frac * lineRect.width);

            // position Now marker
            timelineNow.style.left = `${leftPx}px`;
            timelineNow.style.top = `50%`;
            timelineNow.style.transform = 'translate(-50%, -50%)';

            // position each event item along the same axis
            items.forEach(itemObj => {
                if (!itemObj.date) return;
                const itemFrac = span > 0 ? (itemObj.date.getTime() - minDate.getTime()) / span : 0;
                const itemLeft = lineLeft + (Math.max(0, Math.min(1, itemFrac)) * lineRect.width);
                itemObj.el.style.left = `${itemLeft}px`;
                itemObj.el.style.top = `50%`;
                itemObj.el.style.transform = 'translate(-50%, -50%)';
            });
        } else {
            // vertical: ensure timeline has enough height so we can position items by date
            const minHeight = Math.max(320, Math.round(window.innerHeight * 0.6));
            timeline.style.minHeight = `${minHeight}px`;

            // Force reflow to ensure minHeight is applied before measuring
            void timeline.offsetHeight;

            // recompute rects after ensuring height
            const lineRectV = timelineLine.getBoundingClientRect();
            const timelineRectV = timeline.getBoundingClientRect();

            const lineTop = lineRectV.top - timelineRectV.top;
            const lineHeight = lineRectV.height;
            
            // Position the "Now" marker along the vertical line
            const nowFrac = span > 0 ? (now.getTime() - minDate.getTime()) / span : 0;
            const nowTop = lineTop + (Math.max(0, Math.min(1, nowFrac)) * lineHeight);
            timelineNow.style.position = 'absolute';
            timelineNow.style.left = '50%';
            timelineNow.style.top = `${nowTop}px`;
            timelineNow.style.transform = 'translate(-50%, -50%)';
            
            // Position each event absolutely along the vertical line, alternating left/right
            items.forEach((itemObj, index) => {
                if (!itemObj.date) return;
                const itemFrac = span > 0 ? (itemObj.date.getTime() - minDate.getTime()) / span : 0;
                const itemTop = lineTop + (Math.max(0, Math.min(1, itemFrac)) * lineHeight);
                itemObj.el.style.position = 'absolute';
                itemObj.el.style.top = `${itemTop}px`;
                itemObj.el.style.transform = 'translateY(-50%)';
                
                // alternate left/right by index
                const side = index % 2 === 0 ? 'left' : 'right';
                itemObj.el.setAttribute('data-side', side);
                itemObj.el.style.left = '';
                itemObj.el.style.right = '';
            });
        }
    }

    // initial position and reflow listeners
    positionNowMarker();
    window.addEventListener('resize', () => positionNowMarker());
    // reposition whenever the clock updates (every 30s) in case span changes across midnight
    setInterval(positionNowMarker, 30000);
});
