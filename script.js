// Mobile menu toggle
const menuButton = document.getElementById('menuButton');
const navLinks = document.getElementById('navLinks');

if (menuButton && navLinks) {
    // Toggle menu on click
    menuButton.addEventListener('click', () => {
        menuButton.classList.toggle('active');
        navLinks.classList.toggle('active');
    });

    // Close menu when a link is clicked
    navLinks.querySelectorAll('a').forEach(link => {
        link.addEventListener('click', () => {
            menuButton.classList.remove('active');
            navLinks.classList.remove('active');
        });
    });

    // Close menu when clicking outside
    document.addEventListener('click', (e) => {
        if (!e.target.closest('.navbar')) {
            menuButton.classList.remove('active');
            navLinks.classList.remove('active');
        }
    });
}

// Smooth scroll behavior for navigation links
document.querySelectorAll('a[href^="#"]').forEach(anchor => {
    anchor.addEventListener('click', function (e) {
        e.preventDefault();
        const href = this.getAttribute('href');
        
        // For home link, scroll to top
        if (href === '#home') {
            window.scrollTo({
                top: 0,
                behavior: 'smooth'
            });
            return;
        }
        
        const target = document.querySelector(href);
        if (target) {
            target.scrollIntoView({
                behavior: 'smooth',
                block: 'start'
            });
        }
    });
});

// Hero text fade out on scroll
const heroContent = document.querySelector('.hero-content');
window.addEventListener('scroll', () => {
    const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
    const heroHeight = window.innerHeight;
    
    // Fade out hero content as user scrolls past hero section
    if (scrollTop > 0) {
        const fadeStart = heroHeight * 0.6; // Start fading at 60% of hero height
        const fadeEnd = heroHeight * 1; // Fully fade at end of hero height
        
        if (scrollTop < fadeStart) {
            heroContent.style.opacity = '1';
        } else if (scrollTop > fadeEnd) {
            heroContent.style.opacity = '0';
        } else {
            // Calculate opacity between fadeStart and fadeEnd
            const fadeProgress = (scrollTop - fadeStart) / (fadeEnd - fadeStart);
            heroContent.style.opacity = 1 - fadeProgress;
        }
    }
});

// Add animation on scroll
const observerOptions = {
    threshold: 0.1,
    rootMargin: '0px 0px -100px 0px'
};

const observer = new IntersectionObserver(function(entries) {
    entries.forEach(entry => {
        if (entry.isIntersecting) {
            entry.target.style.opacity = '1';
            entry.target.style.transform = 'translateY(0)';
        }
    });
}, observerOptions);

// Observe cards, items, and section titles for animation
document.querySelectorAll('.about-card, .research-card, .life-card, .cv-item, .research h2, .life h2').forEach(el => {
    el.style.opacity = '0';
    el.style.transform = 'translateY(20px)';
    el.style.transition = 'opacity 0.6s ease, transform 0.6s ease';
    observer.observe(el);
});

// Sticky navbar scroll effect
let lastScrollTop = 0;
const navbar = document.querySelector('.navbar');

window.addEventListener('scroll', () => {
    const scrollTop = window.pageYOffset || document.documentElement.scrollTop;
    
    if (scrollTop > 100) {
        navbar.style.boxShadow = '0 8px 20px rgba(0, 0, 0, 0.2)';
    } else {
        navbar.style.boxShadow = '0 4px 15px rgba(0, 0, 0, 0.15)';
    }
    
    lastScrollTop = scrollTop <= 0 ? 0 : scrollTop;
});

// Highlight active nav link based on scroll position
window.addEventListener('scroll', () => {
    let current = '';
    const sections = document.querySelectorAll('section');
    
    sections.forEach(section => {
        const sectionTop = section.offsetTop;
        const sectionHeight = section.clientHeight;
        if (scrollY >= (sectionTop - 200)) {
            current = section.getAttribute('id');
        }
    });
    
    document.querySelectorAll('.nav-links a').forEach(link => {
        link.style.color = '';
        if (link.getAttribute('href').slice(1) === current) {
            link.style.color = '#f39c12';
        }
    });
});

// Add a subtle parallax effect to hero section
const heroImage = document.querySelector('.hero-image');
if (heroImage) {
    window.addEventListener('scroll', () => {
        const scrollY = window.pageYOffset;
        heroImage.style.transform = `translateY(${scrollY * 0.5}px)`;
    });
}

// Limit scrolling to show only up to the contact section
window.addEventListener('scroll', (e) => {
    const contactSection = document.getElementById('contact');
    if (contactSection) {
        const contactRect = contactSection.getBoundingClientRect();
        const windowHeight = window.innerHeight;
        
        // If contact section is fully visible in viewport, prevent further scrolling
        if (contactRect.bottom <= windowHeight) {
            window.scrollTo(0, contactSection.offsetTop + contactSection.clientHeight - windowHeight);
        }
    }
}, { passive: false });

console.log('Welcome to Richard McManus\'s website!');
