import { ApplicationConfig, provideZoneChangeDetection } from '@angular/core';
import { provideRouter } from '@angular/router';
import {NgcCookieConsentConfig, provideNgcCookieConsent} from 'ngx-cookieconsent';

import { routes } from './app.routes';
import { provideClientHydration, withEventReplay } from '@angular/platform-browser';
import { provideAnimationsAsync } from '@angular/platform-browser/animations/async';
import {provideHttpClient, withFetch} from "@angular/common/http";


export const cookieConfig:NgcCookieConsentConfig = {
  cookie: {
    // TODO: change to proper domain name once we have it
    domain: 'localhost' // or 'your.domain.com' // it is mandatory to set a domain, for cookies to work properly (see https://goo.gl/S2Hy2A)
  },
  palette: {
    popup: {
      background: '#000'
    },
    button: {
      background: '#f1d600'
    }
  },
  theme: 'edgeless',
  type: 'opt-out',
  layout: 'my-custom-layout',
  layouts: {
    "my-custom-layout": '{{messagelink}}{{compliance}}'
  },
  elements:{
    messagelink: `
    <span id="cookieconsent:desc" class="cc-message">{{message}}
      <a aria-label="learn more about our privacy policy" tabindex="1" class="cc-link" href="{{privacyPolicyHref}}"
      target="_blank" rel="noopener">{{privacyPolicyLink}}</a> and our
      <a aria-label="learn more about our terms of service" tabindex="2" class="cc-link" href="{{tosHref}}"
      target="_blank" rel="noopener">{{tosLink}}</a>
    </span>
    `,
  },
  content: {
    message: "This website requires cookies, and the limited processing of your personal data in order to function. " +
      "By using the site you are agreeing to this as outlined in our",

    privacyPolicyLink: 'Privacy Notice',
    privacyPolicyHref: 'https://www.ebi.ac.uk/data-protection/privacy-notice/perturbation-catalogue',

    tosLink: 'Terms of Use',
    tosHref: 'https://www.ebi.ac.uk/about/terms-of-use',
  }
};

export const appConfig: ApplicationConfig = {
  providers: [
    provideZoneChangeDetection({ eventCoalescing: true }),
    provideRouter(routes),
    provideClientHydration(withEventReplay()),
    provideAnimationsAsync(),
    provideHttpClient(
      withFetch()
    ),
    provideNgcCookieConsent(cookieConfig)
  ]
};


