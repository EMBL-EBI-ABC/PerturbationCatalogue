# Perturbation Catalogue front-end

This project was originally generated using [Angular CLI](https://github.com/angular/angular-cli) version 19.0.5.

## Google Cloud Run deployment
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

## Development server

To start a local development server, run:

```bash
ng serve
```

Once the server is running, open your browser and navigate to `http://localhost:4200/`. The application will automatically reload whenever you modify any of the source files.

## Code scaffolding

Angular CLI includes powerful code scaffolding tools. To generate a new component, run:

```bash
ng generate component component-name
```

For a complete list of available schematics (such as `components`, `directives`, or `pipes`), run:

```bash
ng generate --help
```

## Building

To build the project run:

```bash
ng build
```

This will compile your project and store the build artifacts in the `dist/` directory. By default, the production build optimizes your application for performance and speed.

## Running unit tests

To execute unit tests with the [Karma](https://karma-runner.github.io) test runner, use the following command:

```bash
ng test
```

## Running end-to-end tests

For end-to-end (e2e) testing, run:

```bash
ng e2e
```

Angular CLI does not come with an end-to-end testing framework by default. You can choose one that suits your needs.

## Additional Resources

For more information on using the Angular CLI, including detailed command references, visit the [Angular CLI Overview and Command Reference](https://angular.dev/tools/cli) page.
