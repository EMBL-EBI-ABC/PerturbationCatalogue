# Google Cloud Run deployment

1. Go to https://console.cloud.google.com/run.
1. Deploy container → Service → Continuously deploy from a repository (source or function).
1. Set up cloud build.
1. Choose this repository → Next.
1. Build type: Dockerfile; Source location: `/fe/Dockerfile` → Save.
1. Service name: `perturbation-catalogue-fe`.
1. Choose region.
1. Pick: Allow unauthenticated invokations.
1. Pick: CPU is only allocated during request processing.
1. Container(s), volumes, networking, security → Container(s) → Container port: 80
1. Container(s), volumes, networking, security → Container(s) → Variables and Secrets → fill in environment variables:
   - FASTAPI_URL
1. Click: Create.

The deployment can then be accessed at the URL shown on the build page.

Set up path trigger:

1. Go to https://console.cloud.google.com/cloud-build/triggers.
1. Edit the perturbation-catalogue-fe trigger.
1. Click on “Show included and ignored files filters”.
1. Set “Included files filters (glob)” to `fe/**`.
1. Clik on “Save”.

# Auto-generated Angular reference

This project was generated with [Angular CLI](https://github.com/angular/angular-cli) version 18.2.5.

## Development server

Run `ng serve` for a dev server. Navigate to `http://localhost:4200/`. The application will automatically reload if you change any of the source files.

## Code scaffolding

Run `ng generate component component-name` to generate a new component. You can also use `ng generate directive|pipe|service|class|guard|interface|enum|module`.

## Build

Run `ng build` to build the project. The build artifacts will be stored in the `dist/` directory.

## Running unit tests

Run `ng test` to execute the unit tests via [Karma](https://karma-runner.github.io).

## Running end-to-end tests

Run `ng e2e` to execute the end-to-end tests via a platform of your choice. To use this command, you need to first add a package that implements end-to-end testing capabilities.

## Further help

To get more help on the Angular CLI use `ng help` or go check out the [Angular CLI Overview and Command Reference](https://angular.dev/tools/cli) page.
